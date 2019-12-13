from lib.irf_aeff import EffectiveArea
from lib.photometry import Photometrics
from lib.utils import li_ma
import argparse
import copy
import logging
import numpy as np
import os
logger = logging.getLogger('photometry')
# logger.setLevel('DEBUG')
logger.addHandler(logging.StreamHandler())

SOURCE = { 'name': 'Crab', 'ra': 83.6331, 'dec': 22.0145,
           'caldb': 'prod3b-v2', 'irf': 'South_z20_0.5h', }
ENERGY = { 'min': 0.025, 'max': 150.0 }

def counting(filename, pnt, source, rad):
    phm = Photometrics({ 'events_filename': filename })

    off_regions = phm.wobble_regions(pnt, source, rad)
    logging.debug('off regions:'+str(off_regions))

    source_w_rad = copy.deepcopy(source)
    source_w_rad['rad'] = rad
    logging.debug('source obj with rad:'+str(source_w_rad))

    # save regions files in the same dir of events file
    if True:
        tmp_dir, tmp_file = os.path.split(filename)
        phm.write_region([source_w_rad], os.path.join(tmp_dir, 'on.reg'),  color='green', dash=True, width=2)
        phm.write_region(off_regions,    os.path.join(tmp_dir, 'off.reg'), color='red',   dash=True, width=2)

    # counting!!
    on_count = phm.region_counter(source, rad)
    off_count = 0
    for r in off_regions:
        off_count += phm.region_counter(r, r['rad'])
    alpha = 1/len(off_regions)
    excess = on_count - alpha * off_count
    signif = li_ma(on_count, off_count, alpha)
    return on_count, off_count, alpha, excess, signif

# aeff evaluation for source area [cm2]
def eval_aeff(irf_filename, pnt, source, rad, energies, pixel_size):
    aeff = EffectiveArea(irf_filename=irf_filename)
    region = copy.deepcopy(source)
    region['rad'] = rad
    aeff_val = aeff.weighted_value_for_region(region, pnt, energies, pixel_size)
    aeff_val *= 1e4 # cm2
    return aeff_val

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Photometrics analyzer")
    parser.add_argument("events_file", help="the events file")
    parser.add_argument("-sra",  "--source-ra",  help="source right ascension in degrees", type=float, default=None)
    parser.add_argument("-sdec", "--source-dec", help="source declination in degrees",     type=float, default=None)
    parser.add_argument("-pra",  "--pointing-ra",  help="pointing right ascension in degrees", type=float, default=None)
    parser.add_argument("-pdec", "--pointing-dec", help="pointing declination in degrees",     type=float, default=None)
    parser.add_argument("-rad", "--region-radius", help="region radius in degrees", type=float, default=0.2)
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("-irf", "--irf-file", help="the irf file", default=None)
    parser.add_argument("-emin", "--energy-min", help="the energy min to eval the Aeff", type=float, default=None)
    parser.add_argument("-emax", "--energy-max", help="the energy max to eval the Aeff", type=float, default=None)
    parser.add_argument("-psize", "--pixel-size", help="the pixel size to count the Aeff", type=float, default=0.05)
    parser.add_argument("-time", "--livetime", help="time in seconds. Useful only to get the flux [ph/cm²/sec]", type=float, default=None)
    args = parser.parse_args()

    for attr in ['pointing_ra', 'pointing_dec', 'source_ra', 'source_dec']:
        value = getattr(args, attr)
        if value is None:
            raise Exception('The "{}" param cannot be None'.format(attr))

    pnt_coords    = { 'ra': args.pointing_ra, 'dec': args.pointing_dec }
    source_coords = { 'ra': args.source_ra,   'dec': args.source_dec }
    logging.debug('pnt:'+str(pnt_coords)+'src:'+str(source_coords))

    on, off, alpha, excess, significance = counting(args.events_file, pnt_coords, source_coords, args.region_radius)

    source_reg_aeff = eval_aeff(args.irf_file, pnt_coords, source_coords, args.region_radius, [args.energy_min, args.energy_max], args.pixel_size)

    results = [on, off, alpha, excess, significance, source_reg_aeff]
    if args.livetime is not None:
        flux = excess / source_reg_aeff / args.livetime
        results.append(flux)
    print('\t'.join(map(str, results)))

