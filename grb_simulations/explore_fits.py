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

    def simulation_run(self, output_events_file, model_file, time=[0, 1800], energy=[None, None], seed=None, log_file='ctobssim.log', force=False):
        if force or not os.path.isfile(output_events_file):
            sim = ctools.ctobssim()
            sim.clear()
            sim["outevents"] = output_events_file
            sim["inmodel"]   = model_file
            sim["seed"]  = seed if seed else self.seed
            sim["ra"]    = self.ra
            sim["dec"]   = self.dec
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

OBJ = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841, }
ENERGY = { "min": 0.03, "max": 150.0 }

def create_events_file(id, model_file, tmin, tmax):
    working_dir = os.path.dirname(model_file)
    events_file = os.path.join(working_dir, "events_{0:02d}.fits".format(id))
    log_file = os.path.join(working_dir, "ctobssim_{0:02d}.log".format(id))
    obj = Support({ 'name': OBJ['name'], 'ra': OBJ['ra'], 'dec': OBJ['dec'], 'energy_min': ENERGY['min'], 'energy_max': ENERGY['max'], 'seed': None })
    obj.simulation_run(events_file, model_file, time=[tmin, tmax], log_file=log_file)
    return events_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an obeservations list from a spectra fits")
    parser.add_argument("input_fits", help="the fits file with times/energies/spectra data")
    parser.add_argument("-m", "--model", help="the xml template")
    parser.add_argument("-d", "--dir", help="the savings directory (default: data/)", default="data")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("--tmax", help="the final observations time in seconds", type=int)
    args = parser.parse_args()

    # export the model for each time slice
    exporter = TimeSliceExporter(args.input_fits, model_template=args.model, savings_dir=args.dir, verbosity=args.verbose, tmax=args.tmax)
    slices = exporter.export()

    # create events file for each valid time slice
    tstart = 0
    for s in slices:
        tstop = s['tsec']
        if args.tmax is not None:
            if tstart >= args.tmax:
                print("Stop slices loop at slice {}".format(s['slice']), file=sys.stderr)
                break
            if tstop > args.tmax:
                tstop = args.tmax
        s["events_file"] = create_events_file(id=s['slice'], model_file=s['model_file'], tmin=tstart, tmax=tstop)
        tstart=tstop

    # create observation list file
    # http://cta.irap.omp.eu/ctools/users/user_manual/observations.html
    # http://cta.irap.omp.eu/ctools/users/glossary.html#glossary-obsdef
    # http://cta.irap.omp.eu/gammalib/doxygen/classGCTAObservation.html
    obs_list = gammalib.GObservations()
    for s in slices:
        basename, ext = os.path.splitext(os.path.basename(s["model_file"]))
        obs = gammalib.GCTAObservation(s["events_file"])
        obs.id(str(s["slice"]))
        obs.name(basename)
        obs_list.append(obs)
    obs_list_filename = 'obs_definition_list.xml'
    obs_list.save(obs_list_filename)
    exit(0)
