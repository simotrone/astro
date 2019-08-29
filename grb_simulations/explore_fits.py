import argparse
from exporter.timeslice import TimeSliceExporter
import ctools
import gammalib
import os
import sys
import cscripts

class SourceObs:
    def __init__(self, args, verbosity=0):
        self.verbosity = verbosity
        fields = [
            ['name', 'MANDATORY'],
            ['ra',   'MANDATORY'], # source ra
            ['dec',  'MANDATORY'], # source dec
            ['seed',       1],
            ['rad',        5.0],
            ['energy_min', 0.03],
            ['energy_max', 150.0],
            ['caldb',    'prod2'],
            ['irf', 'South_0.5h'],
        ]
        for f in fields:
            field, default = f
            if default == 'MANDATORY':
                self.__dict__[field] = args[field]
            else:
                self.__dict__[field] = args[field] if field in args and args[field] else default

    def simulation_run(self, model_file, events_file, ra=None, dec=None, time=[0, 1800], energy=[None, None], log_file='ctobssim.log', force=False, save=False):
        sim = ctools.ctobssim()
        sim["inmodel"]   = model_file
        sim["outevents"] = events_file
        sim["seed"]  = self.seed
        sim["ra"]    = ra  if ra  else self.ra
        sim["dec"]   = dec if dec else self.dec
        sim["rad"]   = self.rad
        sim["tmin"]  = float(time[0])
        sim["tmax"]  = float(time[1])
        sim["emin"]  = energy[0] if energy[0] else self.energy_min
        sim["emax"]  = energy[1] if energy[1] else self.energy_max
        sim["caldb"] = self.caldb
        sim["irf"]   = self.irf
        sim["logfile"] = log_file
        if force or not os.path.isfile(events_file):
            sim.run()
            sim.logFileOpen()
        else:
            container = gammalib.GObservations()
            gcta_obs = gammalib.GCTAObservation(events_file)
            container.append(gcta_obs)
            sim.obs(container)
        saved = False
        if (save and force) or (save and not os.path.isfile(events_file)):
            sim.save()
            saved = True
        if saved and self.verbosity > 1:
            print("Events file '{}' saved. time [{}-{}]".format(sim["outevents"].value(), tstart, tstop), file=sys.stderr)
        return sim

    def csphagen_run(self, input_obs_list, input_model, source_rad=0.2, output_obs_list='onoff_obs.xml', output_model='onoff_result.xml', log_file='csphagen.log', prefix='onoff', force=False, save=False):
        phagen = cscripts.csphagen()
        if isinstance(input_obs_list, gammalib.GObservations):
            phagen.obs(input_obs_list)
        elif os.path.isfile(input_obs_list):
            # observations list from file
            phagen["inobs"] = input_obs_list
        else:
            raise Exception('Cannot understand input obs list for csphagen')
        phagen["inmodel"] = input_model
        phagen["srcname"] = self.name
        phagen["caldb"]   = self.caldb
        phagen["irf"]     = self.irf
        phagen["ebinalg"]  = "LIN"
        phagen["emin"]     = self.energy_min
        phagen["emax"]     = self.energy_max
        phagen["enumbins"] = 1
        phagen["coordsys"] = "CEL"
        phagen["ra"]    = self.ra
        phagen["dec"]   = self.dec
        phagen["rad"]   = source_rad
        phagen["stack"] = True
        phagen["bkgmethod"] = "REFLECTED"
        phagen["outobs"]   = output_obs_list
        phagen["outmodel"] = output_model
        phagen["prefix"]   = prefix
        phagen["logfile"]  = log_file
        if force or not os.path.isfile(output_obs_list) or not os.path.isfile(output_model):
            phagen.run()
            phagen.logFileOpen()
        elif os.path.isfile(output_obs_list) and os.path.isfile(output_model):
            onoff_obs = gammalib.GObservations(output_obs_list)
            onoff_obs.models(gammalib.GModels(output_model))
            phagen = cscripts.csphagen(onoff_obs)
        else:
            raise Exception("Cannot proceed with csphagen")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_obs_list)) or (save and not os.path.isfile(output_model)):
            phagen.save()
            saved = True
        if saved and self.verbosity > 1:
            print("Files {}, {} created.".format(output_obs_list, output_model), file=sys.stderr)
        return phagen

    def ctlike_run(self, input_obs_list, input_models=None, output_models='ml_result.xml', logfile='ctlike.log', force=False, save=False):
        like = ctools.ctlike()
        if isinstance(input_obs_list, gammalib.GObservations):
            like.obs(input_obs_list)
        elif os.path.isfile(input_obs_list) and os.path.isfile(input_models):
            # observations list from file
            like["inobs"] = input_obs_list
            like["inmodel"]  = input_models
        else:
            raise Exception('Cannot understand input obs list for ctlike')
        like["outmodel"] = result_file
        like["logfile"]  = log_file
        if force or not os.path.isfile(output_models):
            like.run()
            like.logFileOpen()
        elif os.path.isfile(output_models):
            ml_models = gammalib.GModels(output_models)
        else:
            raise Exception("Cannot proceed with ctlike")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_models)):
            like.save()
            saved = True
        if saved and self.verbosity > 1:
            print("File {} created.".format(output_models), file=sys.stderr)
        return like

SOURCE = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841, }
ENERGY = { "min": 0.03, "max": 150.0 }
TIME_SELECTION_SLOTS = [600, 100, 60, 30, 10, 5]

# python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800
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
    sobs = SourceObs({ 'name': SOURCE['name'], 'ra': SOURCE['ra'], 'dec': SOURCE['dec'], 'energy_min': ENERGY['min'], 'energy_max': ENERGY['max'], 'seed': args.seed }, verbosity=args.verbose)

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
    like = sobs.ctlike_run(phagen_obs_list, phagen_obs_list.models(), output_models=like_models_file, logfile=like_log_file, force=args.force)
    if args.verbose > 0:
        print("Maximum Likelihood:\n", like)

    exit(0)
