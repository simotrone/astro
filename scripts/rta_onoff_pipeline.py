import argparse
from lib.utils import li_ma
from lib.photometry import Photometrics
from lib.irf import EffectiveArea
import numpy as np
import math

# Example:
# python rta_onoff_pipeline.py -v -irf test_00_crab/irf_prod3b_v2_South_z20_0.5h.fits -events test_00_crab/events.fits -src-ra 83.6331 -src-dec 22.0145 -pnt-ra 84.1331 -pnt-dec 22.0145 -rad 0.2 -bkgmethod cross --save-off-regions test_00_crab/reflection_off.reg --livetime 1200 -emin 0.025 -emax 150.0 --power-law-index -2.48

def counting(phm, src, rad, off_regions, e_min=None, e_max=None, t_min=None, t_max=None, draconian=False):
    on_count = phm.region_counter(src, rad, emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)
    off_count = 0
    for r in off_regions:
        off_count += phm.region_counter(r, r['rad'], emin=e_min, emax=e_max, tmin=t_min, tmax=t_max)

    alpha = 1 / len(off_regions)
    excess = on_count - alpha * off_count
    signif = li_ma(on_count, off_count, alpha)

    # !!! here we can implement checks
    err_note = None
    if on_count < 10 or off_count < 10:
        err_note = 'Not compliant with Li & Ma requirements.'
        if draconian:
            raise Exception(err_note)

    return on_count, off_count, alpha, excess, signif, err_note

def find_off_regions(phm, algo, src, pnt, rad, verbose=False, save=None):
    if not (src['ra'] and src['dec'] and pnt['ra'] and pnt['dec'] and src['rad']):
        raise Exception('need source and pointing coordinates and a region radius to do aperture photometry')

    off_regions = None
    if algo == 'cross':
        off_regions = phm.cross_regions(pnt, src, rad)
    elif algo == 'reflection':
        off_regions = phm.reflected_regions(pnt, src, rad)
    else:
        raise Exception('invalid background regions algorithm')

    if verbose > 1:
        print('off regions algorithm:', algo)
        for i, o in enumerate(off_regions):
            print('      off regions #{:02d}:'.format(i), o)

    if save:
        phm.write_region(off_regions, save, color='red', dash=True, width=2)

    return off_regions

def aeff_eval(args, src, pnt):
    if not(args.energy_min and args.energy_max and args.pixel_size and args.power_law_index):
        raise Exception('need energy min and max, a pixel size to eval the flux')

    aeff = EffectiveArea(irf_filename=args.irf_file)
    # these IRFs return value in m², so we need convert
    # the source data struct need a 'rad'
    source_reg_aeff = aeff.weighted_value_for_region(src, pnt, [args.energy_min, args.energy_max], args.pixel_size, args.power_law_index) * 1e4 # cm2
    return source_reg_aeff

def main_simple(opts):
    phm = Photometrics({ 'events_filename': opts.events_file })

    src = { 'ra': opts.source_ra,   'dec': opts.source_dec,  'rad': opts.region_radius }
    pnt = { 'ra': opts.pointing_ra, 'dec': opts.pointing_dec }
    radius = opts.region_radius
    livetime = opts.end_time - opts.begin_time

    off_regions = find_off_regions(phm, opts.background_method, src, pnt, radius, verbose=opts.verbose, save=opts.save_off_regions)

    # counting
    on_count, off_count, alpha, excess, significance, err_note = counting(phm, src, radius, off_regions, e_min=opts.energy_min, e_max=opts.energy_max, t_min=opts.begin_time, t_max=opts.end_time, draconian=True)

    flux = None
    if opts.power_law_index:
        region_eff_resp = aeff_eval(opts, src, pnt)
        flux = excess / region_eff_resp / livetime

    ########
    # output
    header=False
    if opts.verbose:
        fmt = '{ra:8.4f} {dec:8.4f} {rad:5.2f} {on:10.1f} {off:10.1f} {alpha:6.3f} {exc:10.2f} {sign:7.2f}'
        output_string = fmt.format(ra=src['ra'], dec=src['dec'], rad=src['rad'], on=on_count, off=off_count, alpha=alpha, exc=excess, sign=significance)

        if not header:
            header_fmt = '{:>8s} {:>8s} {:>5s} {:>10s} {:>10s} {:>6s} {:>10s} {:>7s}'
            header_string = header_fmt.format('src-ra', 'src-dec', 'rad', 'on', 'off', 'alpha', 'excess', 'S')

        if flux:
            header_string += ' '+'{:>15s} {:>17s}'.format('flux ph/cm²/s', 'eff.resp cm²')
            flux_fmt = '{flux:15.3e} {aeff:17.3e}'
            output_string += ' '+flux_fmt.format(flux=flux, aeff=region_eff_resp)

        if not header:
            print(header_string)
            header=True

        print(output_string)

def main_with_time_intervals(opts):
    phm = Photometrics({ 'events_filename': opts.events_file })

    src = { 'ra': opts.source_ra,   'dec': opts.source_dec,  'rad': opts.region_radius }
    pnt = { 'ra': opts.pointing_ra, 'dec': opts.pointing_dec }
    radius = opts.region_radius

    off_regions = find_off_regions(phm, opts.background_method, src, pnt, radius, verbose=opts.verbose, save=opts.save_off_regions)

    # useful to compute the flux. slow exec
    flux = float('NaN')
    region_eff_resp = float('NaN')
    if opts.power_law_index:
        region_eff_resp = aeff_eval(opts, src, pnt)

    # results go in output array
    output = []

    # counter helper
    def counter_fn(t_begin, t_end):
            on_count, off_count, alpha, excess, significance, err_note = counting(phm, src, radius, off_regions, e_min=opts.energy_min, e_max=opts.energy_max, t_min=t_begin, t_max=t_end, draconian=False)

            livetime = t_end - t_begin
            if not math.isnan(region_eff_resp):
                flux = excess / region_eff_resp / livetime

            output.append({ 'on': on_count, 'off': off_count, 'alpha': alpha, 'exc': excess, 'sign': significance, 'err_note': err_note, 'flux': flux, 'aeff': region_eff_resp, 'tmin': t_begin, 'tmax': t_end, 'livetime': livetime })

    if opts.step_time:
        if opts.step_time < 1:
            raise Exception('Step time need to be > 0 sec')
        # stepped counts
        tstart = opts.begin_time
        while(tstart < opts.end_time):
            tstop = tstart + opts.step_time
            if tstop > opts.end_time:
                tstop = opts.end_time
            counter_fn(tstart, tstop)
            tstart = tstop
    else:
        # full counts
        counter_fn(opts.begin_time, opts.end_time)

    print_results(output, src)

########
# output
def print_results(data, src):
    header=False
    for d in data:
        fmt = '{ra:8.4f} {dec:8.4f} {rad:5.2f} {tmin:4.0f}-{tmax:4.0f} {on:10.1f} {off:10.1f} {alpha:6.3f} {exc:10.2f} {sign:7.2f}'
        output_string = fmt.format(ra=src['ra'], dec=src['dec'], rad=src['rad'], **d)

        if not header:
            header_fmt = '{:>8s} {:>8s} {:>5s} {:>9s} {:>10s} {:>10s} {:>6s} {:>10s} {:>7s}'
            header_string = header_fmt.format('src-ra', 'src-dec', 'rad', 'time sec', 'on', 'off', 'alpha', 'excess', 'S')

        if not math.isnan(d['flux']):
            header_string += ' '+'{:>15s} {:>17s} {:>8s}'.format('flux ph/cm²/s', 'eff.resp cm²', 'livetime')
            flux_fmt = '{flux:15.3e} {aeff:17.3e} {livetime:8.1f}'
            output_string += ' '+flux_fmt.format(**d)

        if not header:
            print(header_string)
            header=True

        print(output_string)

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
    parser.add_argument('-bkgmethod', '--background-method', help='choose background regions algorithm. Currently implemented: cross, reflection. (default: cross)', default='cross')
    parser.add_argument('-save-off', '--save-off-regions', help='save off regions in .reg file')
    # aeff options
    parser.add_argument('-emin', '--energy-min', help='the low energy boundary to eval the aeff', type=float)
    parser.add_argument('-emax', '--energy-max', help='the high energy boundary to eval the aeff', type=float)
    parser.add_argument('-psize', '--pixel-size', help="the pixel size to count the Aeff", type=float, default=0.05)
    parser.add_argument('-tbegin', '--begin-time', help="the observation starting time", type=float, default=0)
    parser.add_argument('-tend',   '--end-time', help="the observation ending duration", type=float)
    parser.add_argument('-tstep',  '--step-time', help="the time interval to split the full observation duration", type=float, default=None)
    parser.add_argument('-index', '--power-law-index', help="power law index for aeff calculation", type=float, default=None)

    args = None
    parser_args = parser.parse_args()
    if parser_args.configuration:
        raise Exception('not still implemented')
        args = parsing_xml(parser_args.configuration)
    else:
        args = parser_args

    # main_simple(args)
    main_with_time_intervals(args)

