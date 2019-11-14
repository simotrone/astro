from lib.ctoolswrapper import CToolsWrapper
from lib.exporter.csv import CSVExporter as csvex
from lib.utils import li_ma
from lib.photometry import Photometrics
import numpy as np
import argparse
import logging
import os
import gammalib

logger = logging.getLogger('application')
# logger.setLevel('DEBUG')
logger.addHandler(logging.StreamHandler())

SOURCE = { 'name': 'Crab', 'ra': 83.6331, 'dec': 22.0145,
           'caldb': 'prod3b-v2', 'irf': 'South_z20_0.5h', }
ENERGY = { 'min': 0.025, 'max': 150.0 }
TIME_SELECTION_SLOTS = [600, 100, 60, 30, 20, 10, 5, 4, 3, 2, 1]

parser = argparse.ArgumentParser(description="Create simulations")
parser.add_argument("simulation_model", help="the xml file to simulate initial events")
parser.add_argument("--seed", help="simulations seed (default: 1)", default=1, type=int)
parser.add_argument("--tmax", help="the final observations time in seconds", type=int, default=1800)
parser.add_argument("-ra",  "--ra-shift",  help="simulations right ascension shift in degrees", type=float, default=0)
parser.add_argument("-dec", "--dec-shift", help="simulations declination shift in degrees",     type=float, default=0)
parser.add_argument("-d", "--dir",   help="the savings directory (default: data/)", default="data")
parser.add_argument("-m", "--model", help="source model (for csphagen or ctlike)")
parser.add_argument("--force",     help="force the overwriting", default=False, action="store_true")
parser.add_argument("--save",      help="save the outputs", default=False, action="store_true")
parser.add_argument("-v", "--verbose", action="count", default=0)
args = parser.parse_args()

# prep
sobs = CToolsWrapper({ 'name': SOURCE['name'],
                       'ra':  SOURCE['ra'],
                       'dec': SOURCE['dec'],
                       'energy_min': ENERGY['min'],
                       'energy_max': ENERGY['max'],
                       'caldb': SOURCE['caldb'],
                       'irf':   SOURCE['irf'],
                       'seed': args.seed,
                       'nthreads': 2, }, 
                       verbosity=args.verbose )

working_dir = os.path.join(args.dir, str(args.seed))
try:
    os.makedirs(working_dir)
except FileExistsError as e:
    logger.info("The data dir {} already exists".format(working_dir))

# simulate events
events_file = os.path.join(working_dir, "events.fits")
log_file = os.path.join(working_dir, "ctobssim.log")

pnt = { 'ra':  SOURCE['ra']  + args.ra_shift,
        'dec': SOURCE['dec'] + args.dec_shift }

sim = sobs.simulation_run(model_file=args.simulation_model, events_file=events_file, ra=pnt['ra'], dec=pnt['dec'], time=[0, args.tmax], log_file=log_file, force=args.force, save=args.save)
logger.debug("Simulation:\n"+str(sim))
if sim.obs().size() != 1:
    raise Exception("None or too many simulated observations")

sim_obs_list = sim.obs()

data_to_analyze = []
data_to_analyze.append({ 'tmax': args.tmax,
                         'obs_list': sim_obs_list.clone(),
                         'dir': working_dir,
                         'ra':  pnt['ra'],
                         'dec': pnt['dec'], })

# selections
for t in TIME_SELECTION_SLOTS:
    if t >= args.tmax:
        logger.warning('Skipping time {} because greater of tmax.'.format(t))
        continue

    sel_working_dir = os.path.join(working_dir, "sel_"+str(t))
    if not os.path.isdir(sel_working_dir):
        os.mkdir(sel_working_dir)
    sel_obs_file = os.path.join(sel_working_dir, "events.fits")
    sel_log_file = os.path.join(sel_working_dir, "ctselect.log")
    select = sobs.selection_run(input_obs_list=sim_obs_list, output_obs_list=sel_obs_file, tmin=0, tmax=t, prefix=os.path.join(sel_working_dir, "selected_"), log_file=sel_log_file, force=args.force, save=args.save)
    logger.debug("Selection:\n"+str(select))

    data_to_analyze.append({ 'tmax': t,
                             'obs_list': select.obs().clone(),
                             'dir': sel_working_dir,
                             'ra':  pnt['ra'],
                             'dec': pnt['dec'], })
    logger.info("Selection {} done.".format(sel_working_dir))

def events_gammalib2rec(obs_list):
    events = obs_list[0].events() # GCTAEventList
    fits = gammalib.GFits()
    events.write(fits) # GFits
    events_bintable = fits.table('EVENTS') # GFitsTable
    events_num = events_bintable.nrows()
    tuples = [ (events_bintable['RA'][i], events_bintable['DEC'][i], events_bintable['ENERGY'][i]) for i in range(events_num) ]
    return np.rec.array(tuples, formats='float,float,float', names='RA,DEC,ENERGY')

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
    onoff_log_file = os.path.join(d["dir"], "csphagen.log")
    onoff_obs_file = os.path.join(d["dir"], "onoff_obs_list.xml")
    onoff_model_file = os.path.join(d["dir"], "onoff_result.xml")
    onoff_prefix = os.path.join(d["dir"], "onoff")
    if not args.model:
        raise Exception("Without model cannot run the csphagen process")
    phagen = sobs.csphagen_run(d["obs_list"], input_model=args.model, source_rad=0.2, ebinalg="LOG", enumbins=10, output_obs_list=onoff_obs_file, output_model=onoff_model_file, log_file=onoff_log_file, prefix=onoff_prefix, force=args.force, save=args.save)
    logger.debug("phagen:\n"+str(phagen))
    phagen_obs_list = phagen.obs()
    if phagen_obs_list.size() == 0:
        logger.error("csphagen doesn't provide an on/off observation list for {}/{}".format(d["tmax"], d["dir"]))
        break
    logger.info("on/off {} done.".format(d["dir"]))
    logger.debug("OnOff list:\n"+str(phagen_obs_list))
    logger.debug("GCTAOnOffObservation:\n"+str(phagen_obs_list[0]))
    logger.debug("GCTAOnOffObservation models:\n"+str(phagen_obs_list.models()))

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
    logger.debug("like:\n"+str(like))
    logger.info("maxlike {} done.".format(d["dir"]))
    logger.debug("Maximum Likelihood:\n"+str(like.opt()))
    logger.debug("like observation:\n"+str(like.obs()))
    logger.debug("like models:\n"+str(like.obs().models()))

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
