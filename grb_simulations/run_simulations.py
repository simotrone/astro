import argparse
import gammalib
import os
import sys
from lib.ctoolswrapper import CToolsWrapper
from lib.exporter.csv import CSVExporter as csvex
from lib.utils import Utils.*

SOURCE = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841, }
ENERGY = { "min": 0.03, "max": 150.0 }
TIME_SELECTION_SLOTS = [600, 100, 60, 30, 10, 5]

# python explore_fits.py timeslices.tsv --tmax 1800 --model source_model.xml --dec-shift 0.5 --dir dec_0.5 --save -v
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create simulations and do analysis from a timeslices model")
    parser.add_argument("timeslices_file", help="the tsv file with times/energies/spectra data")
    parser.add_argument("--seed", help="simulations seed (default: 1)", default=1, type=int)
    parser.add_argument("--tmax", help="the final observations time in seconds", type=int, default=1800)
    parser.add_argument("-ra",  "--ra-shift",  help="simulations right ascension shift in degrees", type=float, default=0)
    parser.add_argument("-dec", "--dec-shift", help="simulations declination shift in degrees",     type=float, default=0)
    parser.add_argument("-d", "--dir",   help="the savings directory (default: data/)", default="data")
    parser.add_argument("-m", "--model", help="csphagen model for on-off analysis")
    parser.add_argument("--force",     help="force the overwriting", default=False, action="store_true")
    parser.add_argument("--save",      help="save the outputs", default=False, action="store_true")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    args = parser.parse_args()

    time_slices = read_timeslices_tsv(args.timeslices_file)
    working_dir = os.path.join(args.dir, str(args.seed))
    try:
        os.makedirs(working_dir)
    except FileExistsError as e:
        if args.verbose > 1:
            print("The data dir {} already exists".format(working_dir), file=sys.stderr)

    source = { 'name': SOURCE['name'],
               'ra': SOURCE['ra'],
               'dec': SOURCE['dec'],
               'energy_min': ENERGY['min'],
               'energy_max': ENERGY['max'],
               'seed': args.seed, }
    sobs = CToolsWrapper(source, verbosity=args.verbose)
    sim_obs_list = gammalib.GObservations()
    tstart = 0
    for s in time_slices:
        tstop = float(s['tsec'])
        if args.tmax is not None:
            if tstart >= args.tmax:
                print("Stop time slices loop at slice {}".format(s['id']), file=sys.stderr)
                break
            if tstop > args.tmax:
                tstop = args.tmax
        index = int(s['id'])
        events_file = os.path.join(working_dir, "events_{0:02d}.fits".format(index))
        log_file    = os.path.join(working_dir, "ctobssim_{0:02d}.log".format(index))
        sim = sobs.simulation_run(s['model_file'], events_file, ra=SOURCE['ra']+args.ra_shift, dec=SOURCE['dec']+args.dec_shift, time=[tstart, tstop], log_file=log_file, force=args.force, save=args.save)
        if sim.obs().size() != 1:
            raise Exception("None or too many simulated observations")
        gcta_obs = None
        if args.save:
            # only the files list
            gcta_obs = gammalib.GCTAObservation(events_file)
        else:
            gcta_obs = sim.obs()[0] # GCTAObservation
        gcta_obs.id("{0:02d}".format(index))
        gcta_obs.name("{0}_{1:02d}".format(SOURCE["name"], index))
        sim_obs_list.append(gcta_obs)
        tstart=tstop

    if args.verbose > 0:
        print("Simulations list:\n", sim_obs_list) # GObservations
    if args.save:
        sim_obs_list.save(os.path.join(working_dir, 'sim_obs_list.xml'))

    data_to_analyze = []
    data_to_analyze.append({ 'tmax': args.tmax, 'obs_list': sim_obs_list.clone(), 'dir': working_dir })

    # selections
    for t in TIME_SELECTION_SLOTS:
        sel_working_dir = os.path.join(working_dir, "sel_"+str(t))
        if not os.path.isdir(sel_working_dir):
            os.mkdir(sel_working_dir)
        sel_log_file = os.path.join(sel_working_dir, "ctselect.log")
        sel_obs_file = os.path.join(sel_working_dir, "sel_obs_list.xml")
        select = sobs.selection_run(sim_obs_list, output_obs_list=sel_obs_file, tmin=0, tmax=t, prefix=os.path.join(sel_working_dir, "selected_"), log_file=sel_log_file, force=args.force, save=args.save)
        data_to_analyze.append({ 'tmax': t, 'obs_list': select.obs().clone(), 'dir': sel_working_dir })

    # on/off analysis
    results = []
    for d in data_to_analyze:
        ### csphagen / onoff analysis
        onoff_log_file = os.path.join(d["dir"], "csphagen.log")
        onoff_obs_file = os.path.join(d["dir"], "onoff_obs_list.xml")
        onoff_model_file = os.path.join(d["dir"], "onoff_result.xml")
        onoff_prefix = os.path.join(d["dir"], "onoff")
        phagen = sobs.csphagen_run(d["obs_list"], input_model=args.model, source_rad=0.2, output_obs_list=onoff_obs_file, output_model=onoff_model_file, log_file=onoff_log_file, prefix=onoff_prefix, force=args.force, save=args.save)
        phagen_obs_list = phagen.obs()
        if phagen_obs_list.size() == 0:
            print("csphagen doesn't provide an on/off observation list for {}/{}".format(d["time"], d["dir"]))
            exit(1)
        if args.verbose > 1:
            print("OnOff list:\n", phagen_obs_list)
            print(phagen_obs_list[0]) # GCTAOnOffObservation
            print(phagen_obs_list.models())

        # maximum likelihood against onoff results
        like_models_file = os.path.join(d["dir"], "ml_result.xml")
        like_log_file    = os.path.join(d["dir"], "ctlike.log")
        like = sobs.ctlike_run(phagen_obs_list, input_models=phagen_obs_list.models(), output_models=like_models_file, log_file=like_log_file, force=args.force, save=args.save)
        if args.verbose > 1:
            print("Maximum Likelihood:\n", like.opt())
            print(like.obs())
            print(like.obs().models())

        # summary
        pha_on  = phagen_obs_list[0].on_spec()
        pha_off = phagen_obs_list[0].off_spec()
        ml_models = like.obs().models()
        results.append({ 'name': args.dir,
                         'seed': args.seed,
                         'tmax': d['tmax'],
                         'ts': ml_models[0].ts(),
                         'on_count': pha_on.counts(),
                         'off_count': pha_off.counts(),
                         'alpha': pha_on.backscal(pha_on.size()-1), # spectrum bins-1
                         'li_ma': li_ma(pha_on.counts(), pha_off.counts(), pha_on.backscal(pha_on.size()-1))
                         })

    csvex.save(os.path.join(working_dir, 'results.tsv'), results, headers=list(results[0].keys()), delimiter="\t")
    exit(0)
