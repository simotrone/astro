import argparse
import collections
import logging
import os
from lib.ctoolswrapper import CToolsWrapper

logging.basicConfig(format='%(asctime)s %(levelname)s:\n%(message)s', level=logging.WARNING)

# PYTHONPATH=.. python test_events_generation.py --model ../crab_simulations/crab.xml --save --dir pippo --tmax 100
# PYTHONPATH=../../ python events_generation.py --model ../../crab_simulations/crab.xml --save --dir crab_sim --seed 3 --tmax 1800 --name Crab --ra 83.6331 --dec 22.0145

parser = argparse.ArgumentParser(description="Create simulation")
parser.add_argument("-m", "--model", help="source model for simulation", required=True)
parser.add_argument("-t", "--tmax", help="the final observations time in seconds", type=int, default=1800)
parser.add_argument("--ra",  help="simulations right ascension degrees", type=float, default=0)
parser.add_argument("--dec", help="simulations declination degrees", type=float, default=0)
parser.add_argument("--name", help="source name (descriptive)", default="test")
parser.add_argument("--dir",   help="the savings directory (default: data/)", default="data")
parser.add_argument("--seed",  help="simulations seed (default: 1)", default=1, type=int)
parser.add_argument("--force", help="force the overwriting", default=False, action="store_true")
parser.add_argument("--save",  help="save the outputs", default=False, action="store_true")
parser.add_argument("-v", "--verbose", action="count", default=0)
args = parser.parse_args()

tw = CToolsWrapper({ 'name': args.name,
                     'ra':  args.ra,
                     'dec': args.dec,
                     'energy_min': 0.03, 
                     'energy_max': 150.0,
                     'seed': args.seed, }, verbosity=args.verbose)

working_dir = os.path.join(args.dir, str(args.seed))
try:
    os.makedirs(working_dir)
except FileExistsError as e:
    logging.warning("The data dir {} already exists".format(working_dir))

output_events_filename = os.path.join(working_dir, 'test_events.fits')
log_filename = os.path.join(working_dir, 'test_events_ctobssim.log')

sim = tw.simulation_run( model_file  = args.model,
                         events_file = output_events_filename,
                         time = [0, args.tmax],
                         log_file = log_filename,
                         force = args.force,
                         save  = args.save )

logging.info(sim) # GApplication simulation
logging.info(sim.obs()) # GObservations

g_obs = sim.obs()[0]
print(g_obs) # GObservation
events = g_obs.events()
# events with this attributes:
#   dir() # GCTAInstDir. CTA instrument direction. measured/reconstructed direction of an event
#   energy()
#   event_id() # event identificator 1..N-1
#   index()    # event index 0..N
#   mc_id() # Monte Carlo id, 1 is event from source, 2 is from background
#   phase() # Always 0.0 (?)
#   time()

coords = []
cnt = collections.Counter()
for ev in events:
    # print(ev.event_id(), ev.index(), ev.phase())
    cnt['mc_'+str(ev.mc_id()) ] += 1
    cnt['total'] += 1
    inst_dir = ev.dir()
    coords.append((inst_dir.theta(), inst_dir.phi()))

print("Events counter:", cnt)

# show details about http://cta.irap.omp.eu/gammalib/doxygen/classGCTAEventAtom.html
def show_event_details(ev):
    # http://cta.irap.omp.eu/gammalib/doxygen/classGCTAInstDir.html
    gcta_inst_dir = ev.dir()
    # http://cta.irap.omp.eu/gammalib/doxygen/classGSkyDir.html
    g_sky_dir = gcta_inst_dir.dir()
    # http://cta.irap.omp.eu/gammalib/doxygen/classGTime.html
    g_time = ev.time()
    print("""
        event: {}
        time: {} sec, {} jd, {} mjd, {} utc
        GCTAInstDir: {}
        GSkyDir (gal deg): {} {}
        GSkyDir (gal rad): {} {}
        GSkyDir (cel deg): {} {}
        GSkyDir (cel rad): {} {}
        """.format(ev,
            g_time, g_time.jd(), g_time.mjd(), g_time.utc(),
            gcta_inst_dir,
            g_sky_dir.l_deg(), g_sky_dir.b_deg(),
            g_sky_dir.l(), g_sky_dir.b(),
            g_sky_dir.ra_deg(), g_sky_dir.dec_deg(),
            g_sky_dir.ra(), g_sky_dir.dec()))
    return

show_event_details(events[0])

