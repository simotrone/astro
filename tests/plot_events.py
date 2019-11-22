import zlib
import argparse
import collections
import gammalib
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

logging.basicConfig(format='%(asctime)s %(levelname)s:\n%(message)s', level=logging.WARNING)

# PYTHONPATH=.. python test_events_generation.py --model ../crab_simulations/crab.xml --save --dir pippo --tmax 100
# PYTHONPATH=../../ python events_generation.py --model ../../crab_simulations/crab.xml --save --dir crab_sim --seed 3 --tmax 1800 --name Crab --ra 83.6331 --dec 22.0145

parser = argparse.ArgumentParser(description="Plot events")
parser.add_argument("input_file", help="the input events file")
parser.add_argument("--dir",   help="the savings directory (default: data/)", default="data")
parser.add_argument("--force", help="force the overwriting", default=False, action="store_true")
parser.add_argument("--save",  help="save the outputs", default=False, action="store_true")
parser.add_argument("-v", "--verbose", action="count", default=0)
args = parser.parse_args()

# events with this attributes:
#   dir() # GCTAInstDir. CTA instrument direction. measured/reconstructed direction of an event
#   energy()
#   event_id() # event identificator 1..N-1
#   index()    # event index 0..N
#   mc_id() # Monte Carlo id, 1 is event from source, 2 is from background
#   phase() # Always 0.0 (?)
#   time()

events = []
try:
    events_list = gammalib.GCTAEventList(args.input_file) # http://cta.irap.omp.eu/gammalib/doxygen/classGCTAEventList.html
    for ev in events_list:
        events.append(ev)
except:
    g_obss = gammalib.GObservations(args.input_file)
    for gcta_obs in g_obss:
        if gcta_obs.has_events():
            events_list = gcta_obs.events() # GCTAEventList
            for ev in events_list:
                events.append(ev)

data = { 'l':[], 'b':[], 'ene':[], 'mc_id':[], 'detx':[], 'dety':[] }
cnt = collections.Counter()
for ev in events:
    gcta_inst_dir = ev.dir()
    g_sky_dir = gcta_inst_dir.dir() # http://cta.irap.omp.eu/gammalib/doxygen/classGSkyDir.html
    data['l'].append(g_sky_dir.l_deg())
    data['b'].append(g_sky_dir.b_deg())

    data['ene'].append(ev.energy().TeV())

    data['mc_id'].append(ev.mc_id())
    cnt['mc_'+str(ev.mc_id())] += 1
    cnt['total'] += 1

    data['detx'].append(gcta_inst_dir.detx())
    data['dety'].append(gcta_inst_dir.dety())

print("Events counter:", cnt)

PLOT = True
######
# Plot N skycoord
if PLOT:
    LIMIT = 1000
    plt.scatter(data['l'][:LIMIT], data['b'][:LIMIT], s=10, alpha=0.5)
    plt.xlabel("Gal long (deg)")
    plt.ylabel("Gal lat (deg)")
    plt.title('First {} events'.format(LIMIT))
    plt.show()

######
# Plot Energy
ene_bins = np.logspace(-2,2, num=100)
if PLOT:
    plt.hist(data['ene'], bins=ene_bins)
    plt.semilogx()
    plt.xlabel("Event energy (TeV)")
    plt.ylabel("Number of events")
    plt.show()

######
# Plot gamma vs hadron (?)
# From here: https://docs.gammapy.org/0.8/notebooks/cta_1dc_introduction.html
# not sure about this: mc_id is events(#1) vs background events(#2)
# Furthermore I do graph
print("Number of events:", cnt['total'])
print("Number of gammas:", cnt['mc_2'])
print("Number of hadrons:", cnt['mc_1'])

if PLOT:
    ene_gamma  = []
    ene_hadron = []
    for i, mc in enumerate(data['mc_id']):
        if mc != 1:
            ene_gamma.append(data['ene'][i])
        else:
            ene_hadron.append(data['ene'][i])
    
    opts = dict(bins=ene_bins, density=True, histtype="step")
    plt.hist(ene_hadron, label='hadron (mc == 1)', **opts)
    plt.hist(ene_gamma, label='gamma (mc != 1)', **opts)
    plt.loglog()
    plt.xlabel('Event energy (TeV)')
    plt.ylabel('Number of events')
    plt.legend()
    plt.show()

if PLOT:
    offset = np.sqrt(np.array(data['detx']) **2 + np.array(data['dety']) **2)
    offset_bins = np.arange(offset.min(), offset.max(), 0.001)
    hist, xedges, yedges = np.histogram2d(x=data['ene'], y=offset, bins=(ene_bins, offset_bins))

    from matplotlib.colors import LogNorm
    plt.pcolormesh(ene_bins, offset_bins, hist.T, norm=LogNorm())
    plt.semilogx()
    plt.colorbar()
    plt.xlabel("Energy (TeV)")
    plt.ylabel("Offset (deg)")
    plt.show()


exit(0)

# show details about http://cta.irap.omp.eu/gammalib/doxygen/classGCTAEventAtom.html

