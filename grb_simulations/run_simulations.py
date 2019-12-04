import argparse
import gammalib
import logging
import os
import sys
from lib.ctoolswrapper import CToolsWrapper
from lib.exporter.csv import CSVExporter as csvex
from lib.utils import li_ma, read_timeslices_tsv
from lib.photometry import Photometrics
import numpy as np

SOURCE = { "name": "run0406_ID000126", "ra": 33.057, "dec": -51.841,
           'caldb': 'prod3b-v2', 'irf': 'South_z40_0.5h', }
ENERGY = { "min": 0.025, "max": 150.0 }
TIME_SELECTION_SLOTS = [600, 100, 60, 30, 20, 10, 5, 4, 3, 2, 1]

# python explore_fits.py timeslices.tsv --tmax 1800 --model source_model.xml --dec-shift 0.5 --dir dec_0.5 --save -v
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
    logging.warning("The data dir {} already exists".format(working_dir))

sobs = CToolsWrapper({ 'name': SOURCE['name'],
                       'ra': SOURCE['ra'],
                       'dec': SOURCE['dec'],
                       'energy_min': ENERGY['min'],
                       'energy_max': ENERGY['max'],
                       'caldb': SOURCE['caldb'],
                       'irf':   SOURCE['irf'],
                       'seed': args.seed,
                       'nthreads': 2, },
                       verbosity=args.verbose )
pnt = { 'ra':  SOURCE['ra']  + args.ra_shift,
        'dec': SOURCE['dec'] + args.dec_shift }

sim_obs_list = gammalib.GObservations()
tstart = 0
for s in time_slices:
    tstop = float(s['tsec'])
    if args.tmax is not None:
        if tstart >= args.tmax:
            logging.info("Stop time slices loop at slice {} where time >= limit".format(s['id']))
            break
        if tstop > args.tmax:
            tstop = args.tmax
    index = int(s['id'])
    events_file = os.path.join(working_dir, "events_{0:02d}.fits".format(index))
    log_file    = os.path.join(working_dir, "ctobssim_{0:02d}.log".format(index))
    sim = sobs.simulation_run(s['model_file'], events_file, ra=pnt['ra'], dec=pnt['dec'], time=[tstart, tstop], log_file=log_file, force=args.force, save=args.save)
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
    logging.warning("Simulation {} done.".format(gcta_obs.name()))
    tstart=tstop

if args.verbose > 0:
    print("Simulations list:\n", sim_obs_list) # GObservations
if args.save:
    sim_obs_list.save(os.path.join(working_dir, 'sim_obs_list.xml'))

data_to_analyze = []
data_to_analyze.append({ 'tmax': args.tmax,
                         'obs_list': sim_obs_list.clone(),
                         'dir': working_dir,
                         'ra':  pnt['ra'],
                         'dec': pnt['dec'],
                         'emin': ENERGY['min'],
                         'emax': ENERGY['max'], })

# selections
for t in TIME_SELECTION_SLOTS:
    if t >= args.tmax:
        logging.warning('Skipping time {} because greater of tmax.'.format(t))
        continue

    sel_working_dir = os.path.join(working_dir, 'sel_{0:04d}_{1:.3f}_{2:.3f}'.format(t, ENERGY['min'], ENERGY['max']))
    if not os.path.isdir(sel_working_dir):
        os.mkdir(sel_working_dir)
    sel_log_file = os.path.join(sel_working_dir, "ctselect.log")
    sel_obs_file = os.path.join(sel_working_dir, "sel_obs_list.xml")
    select = sobs.selection_run(sim_obs_list, output_obs_list=sel_obs_file, tmin=0, tmax=t, prefix=os.path.join(sel_working_dir, "selected_"), log_file=sel_log_file, force=args.force, save=args.save)
    data_to_analyze.append({ 'tmax': t,
                             'obs_list': select.obs().clone(),
                             'dir': sel_working_dir,
                             'ra':  pnt['ra'],
                             'dec': pnt['dec'],
                             'emin': ENERGY['min'],
                             'emax': ENERGY['max'], })
    logging.info("Selection {} done.".format(sel_working_dir))

def events_gammalib2rec(obs_list):
    ev_list = []
    for obs in obs_list:
        events = obs.events() # GCTAEventList
        fits = gammalib.GFits()
        events.write(fits) # GFits
        events_bintable = fits.table('EVENTS') # GFitsTable
        events_num = events_bintable.nrows()
        tuples = [ (events_bintable['RA'][i], events_bintable['DEC'][i], events_bintable['ENERGY'][i]) for i in range(events_num) ]
        ev_list += tuples
    return np.rec.array(ev_list, formats='float,float,float', names='RA,DEC,ENERGY')

def photometrics_counts(data):
    events_list = events_gammalib2rec(data['obs_list'])
    phm = Photometrics({ 'events_list': events_list })
    pnt_coords = { 'ra': data['ra'], 'dec': data['dec'] }
    source_coords = { 'ra': SOURCE['ra'], 'dec': SOURCE['dec'] }
    region_rad = 0.2
    reflected_regions = phm.reflected_regions(pnt_coords, source_coords, region_rad)
    on_count = phm.region_counter(source_coords, region_rad) # no energy thresholds, emin=data['emin'], emax=data['emax'])
    off_count = 0
    for r in reflected_regions:
        off_count += phm.region_counter(r, r['rad']) # no energy thresholds emin=data['emin'], emax=data['emax'])
    alpha = 1/len(reflected_regions)
    return { 'on': on_count, 'off': off_count, 'alpha': alpha, 'excess': on_count - alpha * off_count }

# on/off analysis
results = []
for d in data_to_analyze:
    photometrics_results = photometrics_counts(d)

    ### csphagen / onoff analysis
    onoff_log_file = os.path.join(d['dir'], 'csphagen.log')
    onoff_obs_file = os.path.join(d['dir'], 'onoff_obs_list.xml')
    onoff_model_file = os.path.join(d['dir'], 'onoff_result.xml')
    onoff_prefix = os.path.join(d['dir'], 'onoff')
    if not args.model:
        raise Exception('Without model cannot run the csphagen process')
    phagen = sobs.csphagen_run(d['obs_list'], input_model=args.model, source_rad=0.2, ebinalg="LOG", enumbins=1, output_obs_list=onoff_obs_file, output_model=onoff_model_file, log_file=onoff_log_file, prefix=onoff_prefix, stacked=True, force=args.force, save=args.save)
    phagen_obs_list = phagen.obs()
    if phagen_obs_list.size() == 0:
        logging.error("csphagen doesn't provide an on/off observation list for {}/{}".format(d['tmax'], d['dir']))
        break
    logging.info('on/off {} done.'.format(d['dir']))
    logging.debug('OnOff list:\n', phagen_obs_list)
    logging.debug(phagen_obs_list[0]) # GCTAOnOffObservation
    logging.debug(phagen_obs_list.models())

    pha_on  = phagen_obs_list[0].on_spec()
    pha_off = phagen_obs_list[0].off_spec()
    on_count  = pha_on.counts()
    off_count = pha_off.counts()
    alpha = pha_on.backscal(pha_on.size()-1) # spectrum bins-1
    excess_count = on_count - alpha * off_count

    # maximum likelihood against onoff results
    like_models_file = os.path.join(d["dir"], "ml_result.xml")
    like_log_file    = os.path.join(d["dir"], "ctlike.log")
    like = sobs.ctlike_run(phagen_obs_list, input_models=phagen_obs_list.models(), output_models=like_models_file, log_file=like_log_file, force=args.force, save=args.save)
    logging.info("maxlike {} done.".format(d["dir"]))
    logging.debug("Maximum Likelihood:\n"+str(like.opt()))
    logging.debug("like observation:\n"+str(like.obs()))
    logging.debug("like models:\n"+str(like.obs().models()))

    # summary
    ml_models = like.obs().models()
    if args.force or not os.path.isfile(like_models_file):
        ml_models.save(like_models_file)

    spectral_model = ml_models[0].spectral()
    GENERGY = { 'min': gammalib.GEnergy(ENERGY['min'], 'TeV'),
                'max': gammalib.GEnergy(ENERGY['max'], 'TeV'), }
    flux  = spectral_model.flux(GENERGY['min'], GENERGY['max'])  # ph/cm²/s
    eflux = spectral_model.eflux(GENERGY['min'], GENERGY['max']) # erg/cm²/s

    results.append({ 'name': args.dir,
                     'ra':   pnt['ra'],
                     'dec':  pnt['dec'],
                     'seed': args.seed,
                     'tmax': d['tmax'],
                     'emin': d['emin'],
                     'emax': d['emax'],
                     'ts': ml_models[0].ts(),
                     'index_value': ml_models[0]['Index'].value(),
                     'index_error': ml_models[0]['Index'].error(),
                     'prefactor_value': ml_models[0]['Prefactor'].value(),
                     'prefactor_error': ml_models[0]['Prefactor'].error(),
                     'pivot_value': ml_models[0]['PivotEnergy'].value(),
                     'pivot_error': ml_models[0]['PivotEnergy'].error(),
                     'flux': flux,
                     'eflux': eflux,
                     'on_count': on_count,
                     'off_count': off_count,
                     'alpha': alpha,
                     'excess_count': excess_count,
                     'li_ma': li_ma(on_count, off_count, alpha),
                     'phm_on': photometrics_results['on'],
                     'phm_off': photometrics_results['off'],
                     'phm_excess': photometrics_results['excess'],
                     'phm_li_ma': li_ma(photometrics_results['on'], photometrics_results['off'], photometrics_results['alpha']),
                     })

try:
	csvex.save(os.path.join(args.dir, 'results_{}.tsv'.format(str(args.seed))), results, headers=list(results[0].keys()), delimiter="\t")
except:
	print(results)

exit(0)
