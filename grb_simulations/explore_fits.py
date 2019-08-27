import argparse
from exporter.timeslice import TimeSliceExporter
import ctools
import gammalib
import os
import sys

class Support:
    def __init__(self, args):
        fields = [
            ['name', 'MANDATORY'],
            ['ra',   'MANDATORY'],
            ['dec',  'MANDATORY'],
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

    def simulation_run(self, output_events_file, model_file, ra=None, dec=None, time=[0, 1800], energy=[None, None], log_file='ctobssim.log', force=False):
        if force or not os.path.isfile(output_events_file):
            sim = ctools.ctobssim()
            sim.clear()
            sim["outevents"] = output_events_file
            sim["inmodel"]   = model_file
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
            sim.logFileOpen()
            sim.run()
            sim.save()
            print("Events file '{}' time [{}-{}]".format(output_events_file, time[0], time[1]), file=sys.stderr)
        return output_events_file

# http://cta.irap.omp.eu/ctools/users/user_manual/observations.html
# http://cta.irap.omp.eu/ctools/users/glossary.html#glossary-obsdef
# http://cta.irap.omp.eu/gammalib/doxygen/classGCTAObservation.html
def create_observations_list(ref, output_filename):
    obs_list = gammalib.GObservations()
    for s in ref:
        basename, ext = os.path.splitext(os.path.basename(s["model_file"]))
        obs = gammalib.GCTAObservation(s["events_file"])
        obs.id(str(s["id"]))
        obs.name(basename)
        obs_list.append(obs)
    obs_list.save(output_filename)
    return obs_list

def create_events_selections(observations_list, time_slots):
    return

SOURCE = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841, }
ENERGY = { "min": 0.03, "max": 150.0 }
TIME_SLOTS = [1800, 600, 100, 60, 30, 10, 5]

# python explore_fits.py run0406_ID000126.fits --template-model run0406_ID000126.xml --tmax 1800
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an obeservations list from a spectra fits")
    parser.add_argument("input_fits",          help="the fits file with times/energies/spectra data")
    parser.add_argument("-ra",  "--ra-shift",  help="right ascension shift in degrees", type=float, default=0)
    parser.add_argument("-dec", "--dec-shift", help="declination shift in degrees",     type=float, default=0)
    parser.add_argument("-tm",  "--template-model", help="the xml template")
    parser.add_argument("-d", "--dir", help="the savings directory (default: data/)", default="data")
    parser.add_argument("--tmax",      help="the final observations time in seconds", type=int)
    parser.add_argument("--force",     help="force the overwriting", default=False, action="store_true")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("--seed", default=1, type=int)
    args = parser.parse_args()

    # export the model for each time slice
    exporter = TimeSliceExporter(args.input_fits, template_model=args.template_model, savings_dir=args.dir, verbosity=args.verbose, tmax=args.tmax)
    slices = exporter.export(force=args.force)

    # create events file for each valid time slice
    obj = Support({ 'name': SOURCE['name'], 'ra': SOURCE['ra'], 'dec': SOURCE['dec'], 'energy_min': ENERGY['min'], 'energy_max': ENERGY['max'], 'seed': args.seed })
    tstart = 0
    for s in slices:
        tstop = s['tsec']
        if args.tmax is not None:
            if tstart >= args.tmax:
                print("Stop slices loop at slice {}".format(s['id']), file=sys.stderr)
                break
            if tstop > args.tmax:
                tstop = args.tmax
        working_dir = os.path.dirname(s['model_file'])
        events_file = os.path.join(working_dir, "events_{0:02d}.fits".format(s['id']))
        log_file    = os.path.join(working_dir, "ctobssim_{0:02d}.log".format(s['id']))
        s["events_file"] = obj.simulation_run(events_file, model_file=s['model_file'], ra=SOURCE['ra']+args.ra_shift, dec=SOURCE['dec']+args.dec_shift, time=[tstart, tstop], log_file=log_file, force=args.force)
        tstart=tstop

    obs_list_filename = '{}_obs_definition_list.xml'.format(args.dir)
    obs_list = create_observations_list(ref=slices, output_filename=obs_list_filename)
    if args.verbose > 0:
        print(obs_list)

    # selections
    create_events_selections(obs_list_filename, time_slots=TIME_SLOTS)

    # binning
    # csphagen
    # obj.csphagen_run(obs_list_filename)
    # likelihood
    exit(0)

