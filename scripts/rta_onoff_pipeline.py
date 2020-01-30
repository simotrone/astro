import argparse
from lib.utils import li_ma
from lib.photometry import Photometrics
from lib.irf import EffectiveArea
import numpy as np

# Example:
# python rta_onoff_pipeline.py -v -irf test_00_crab/irf_prod3b_v2_South_z20_0.5h.fits -events test_00_crab/events.fits -src-ra 83.6331 -src-dec 22.0145 -pnt-ra 84.1331 -pnt-dec 22.0145 -rad 0.2 -bkgmethod wobble --save-off-regions test_00_crab/reflection_off.reg --livetime 1200 -emin 0.025 -emax 150.0 --power-law-index -2.48

def counting(phm, src, rad, off_regions, verbose=False):
    on_count = phm.region_counter(src, rad)
    off_count = 0
    for r in off_regions:
        off_count += phm.region_counter(r, r['rad'])

    alpha = 1 / len(off_regions)
    excess = on_count - alpha * off_count
    signif = li_ma(on_count, off_count, alpha)

    if verbose:
        print('counting for src:', src)
        print('        on count:', on_count)
        print('       off count:', off_count)
        print('           alpha: {:.3f}'.format(alpha))
        print('          excess: {:.2f}'.format(excess))
        print('    significance: {:.2f}'.format(signif))
    return on_count, off_count, alpha, excess, signif

def find_off_regions(phm, algo, src, pnt, rad, verbose=False, save=None):
    off_regions = None
    if algo == 'wobble':
        off_regions = phm.wobble_regions(pnt, src, rad)
    elif algo == 'reflection':
        off_regions = phm.reflected_regions(pnt, src, rad)
    else:
        raise Exception('invalid background regions algorithm')

    if verbose:
        print('off regions algorithm:', algo)
        for i, o in enumerate(off_regions):
            print('      off regions #{:02d}:'.format(i), o)

    if save:
        phm.write_region(off_regions, save, color='red', dash=True, width=2)

    return off_regions

def main(opts):
    aeff = EffectiveArea(irf_filename=opts.irf_file)
    phm = Photometrics({ 'events_filename': opts.events_file })

    # regions
    if not (opts.source_ra and opts.source_dec and opts.pointing_ra and opts.pointing_dec and opts.region_radius):
        raise Exception('need source and pointing coordinates and a region radius to do aperture photometry')

    src = { 'ra': opts.source_ra,   'dec': opts.source_dec,  'rad': opts.region_radius }
    pnt = { 'ra': opts.pointing_ra, 'dec': opts.pointing_dec }
    radius = opts.region_radius

    off_regions = find_off_regions(phm, opts.background_method, src, pnt, radius, verbose=opts.verbose, save=opts.save_off_regions)

    # counting
    on_count, off_count, alpha, excess, significance = counting(phm, src, radius, off_regions, verbose=opts.verbose)

    # !!! here we can implement checks
    if on_count < 10 or off_count < 10:
        raise Exception('Not compliant with Li & Ma requirements.')

    # flux
    if not(opts.energy_min and opts.energy_max and opts.pixel_size and opts.livetime and opts.power_law_index):
        raise Exception('need energy min and max, a pixel size and livetime to eval the flux')

    # these IRFs return value in m², so we need convert
    # the source data struct need a 'rad'
    source_reg_aeff = aeff.weighted_value_for_region(src, pnt, [opts.energy_min, opts.energy_max], opts.pixel_size, opts.power_law_index) * 1e4 # cm2
    livetime = opts.livetime
    flux = excess / source_reg_aeff / livetime
    if opts.verbose:
        print('  effective area: {:.3e} cm²'.format(source_reg_aeff))
        print('        livetime:', livetime)
        print('            flux: {:.3e} ph/cm²/s'.format(flux))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="on off analysis from an event list")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument('--configuration', help='read xml configuration for parameters')
    parser.add_argument('-events', '--events-file', help='the events file (.fits)')
    parser.add_argument("-irf", "--irf-file", help="the irf file", default=None)
    # regions options
    parser.add_argument('-src-ra', '--source-ra', help='the source right ascension', type=float)
    parser.add_argument('-src-dec', '--source-dec', help='the source declination', type=float)
    parser.add_argument('-pnt-ra', '--pointing-ra', help='the pointing right ascension', type=float)
    parser.add_argument('-pnt-dec', '--pointing-dec', help='the pointing declination', type=float)
    parser.add_argument('-rad', '--region-radius', help='the region radius (default: 0.2°)', default=0.2, type=float)
    parser.add_argument('-bkgmethod', '--background-method', help='choose background regions algorithm. Currently implemented: wobble, reflection. (default: wobble)', default='wobble')
    parser.add_argument('-save-off', '--save-off-regions', help='save off regions in .reg file')
    # aeff options
    parser.add_argument('-emin', '--energy-min', help='the low energy boundary to eval the aeff', type=float)
    parser.add_argument('-emax', '--energy-max', help='the high energy boundary to eval the aeff', type=float)
    parser.add_argument('-psize', '--pixel-size', help="the pixel size to count the Aeff", type=float, default=0.05)
    parser.add_argument('-t', '--livetime', help="the observation duration", type=float)
    parser.add_argument('-index', '--power-law-index', help="power law index for aeff calculation", type=float, default=-2.40)

    args = None
    parser_args = parser.parse_args()
    if parser_args.configuration:
        raise Exception('not still implemented')
        args = parsing_xml(parser_args.configuration)
    else:
        args = parser_args

    main(args)

