import argparse
import gammalib
import os
import sys
from lib.exporter.timeslice import TimeSliceExporter
from lib.ctoolswrapper import CToolsWrapper

SOURCE = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841, }
ENERGY = { "min": 0.03, "max": 150.0 }
TIME_SELECTION_SLOTS = [600, 100, 60, 30, 10, 5]

# python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800
# python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800 --model source_model.xml --dec-shift 0.5 --dir dec_0.5 --save -v
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an obeservations list from a spectra fits")
    parser.add_argument("input_fits",          help="the fits file with times/energies/spectra data")
    parser.add_argument("-ra",  "--ra-shift",  help="right ascension shift in degrees", type=float, default=0)
    parser.add_argument("-dec", "--dec-shift", help="declination shift in degrees",     type=float, default=0)
    parser.add_argument("-tm",  "--template-model", help="the xml template")
    parser.add_argument("-m",   "--model",     help="csphagen model for on-off analysis")
    parser.add_argument("--tmax",      help="the final observations time in seconds", type=int, default=1800)
    parser.add_argument("-d", "--dir", help="the savings directory (default: data/)", default="data")
    parser.add_argument("--force",     help="force the overwriting", default=False, action="store_true")
    parser.add_argument("--save",      help="save the outputs", default=False, action="store_true")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("--seed", default=1, type=int)
    args = parser.parse_args()

    # export the model for each time slice
    exporter = TimeSliceExporter(args.input_fits, template_model=args.template_model, savings_dir=args.dir, verbosity=args.verbose, tmax=args.tmax)
    time_slices = exporter.export(force=args.force)

    working_dir = os.path.join(args.dir, str(args.seed))
    try:
        os.makedirs(working_dir)
    except FileExistsError as e:
        if args.verbose > 1:
            print("The data dir already exists", file=sys.stderr)

    # create events file for each valid time slice
    sobs = CToolsWrapper({ 'name': SOURCE['name'], 'ra': SOURCE['ra'], 'dec': SOURCE['dec'], 'energy_min': ENERGY['min'], 'energy_max': ENERGY['max'], 'seed': args.seed }, verbosity=args.verbose)

    # http://cta.irap.omp.eu/ctools/users/user_manual/observations.html
    # http://cta.irap.omp.eu/ctools/users/glossary.html#glossary-obsdef
    # http://cta.irap.omp.eu/gammalib/doxygen/classGCTAObservation.html
    sim_obs_list = gammalib.GObservations()
    tstart = 0
    for s in time_slices:
        tstop = s['tsec']
        if args.tmax is not None:
            if tstart >= args.tmax:
                print("Stop time slices loop at slice {}".format(s['id']), file=sys.stderr)
                break
            if tstop > args.tmax:
                tstop = args.tmax
        events_file = os.path.join(working_dir, "events_{0:02d}.fits".format(s['id']))
        log_file    = os.path.join(working_dir, "ctobssim_{0:02d}.log".format(s['id']))
        sim = sobs.simulation_run(s['model_file'], events_file, ra=SOURCE['ra']+args.ra_shift, dec=SOURCE['dec']+args.dec_shift, time=[tstart, tstop], log_file=log_file, force=args.force, save=args.save)
        if sim.obs().size() != 1:
            raise Exception("None or too many simulated observations")
        gcta_obs = sim.obs()[0] # GCTAObservation
        gcta_obs.id("{0:02d}".format(s["id"]))
        gcta_obs.name("{0}_{1:02d}".format(SOURCE["name"], s["id"]))
        sim_obs_list.append(gcta_obs)
        tstart=tstop

    if args.verbose > 0:
        print("Simulations list:\n", sim_obs_list) # GObservations
    if args.save:
        sim_obs_list.save(os.path.join(working_dir, 'sim_obs_list.xml'))

    # selections

    ### csphagen / onoff analysis
    onoff_log_file = os.path.join(working_dir, "csphagen.log")
    onoff_obs_file = os.path.join(working_dir, "onoff_obs_list.xml")
    onoff_model_file = os.path.join(working_dir, "onoff_result.xml")
    prefix = os.path.join(working_dir, "onoff")
    phagen = sobs.csphagen_run(sim_obs_list, input_model=args.model, source_rad=0.2, output_obs_list=onoff_obs_file, output_model=onoff_model_file, log_file=onoff_log_file, prefix=prefix, force=args.force, save=args.save)
    phagen_obs_list = phagen.obs()
    if args.verbose > 0:
        print("OnOff list:\n", phagen_obs_list)
        print(phagen_obs_list[0]) # GCTAOnOffObservation
        print(phagen_obs_list.models())

    # maximum likelihood
    like_models_file = os.path.join(working_dir, "ml_result.xml")
    like_log_file    = os.path.join(working_dir, "ctlike.log")
    like = sobs.ctlike_run(phagen_obs_list, input_models=phagen_obs_list.models(), output_models=like_models_file, log_file=like_log_file, force=args.force, save=args.save)
    # like = sobs.ctlike_run(onoff_obs_file, input_models=onoff_model_file, output_models=like_models_file, log_file=like_log_file, force=args.force)
    if args.verbose > 0:
        print("Maximum Likelihood:\n", like)
        print(like.obs())
        print(like.obs().models())


    exit(0)
