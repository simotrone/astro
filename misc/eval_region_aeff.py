from lib.irf import EffectiveArea
from lib import utils
from astropy.coordinates import SkyCoord
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
# a better estimation of effective area in region

def plot_region(region, points, pointing=None, title=None):
    fig, ax = plt.subplots()

    circle = plt.Circle((region['ra'], region['dec']), region['rad'], color='g', linestyle='--', fill=False)
    ax.add_artist(circle)
    ax.text(region['ra']-region['rad'], region['dec']-region['rad']-0.08, 'Source region')

    midpoints_ra = []
    midpoints_dec = []
    for p in points:
        midpoints_ra.append(p['ra'])
        midpoints_dec.append(p['dec'])
    ax.plot(midpoints_ra, midpoints_dec, 'rx', markersize=2)
    # ax.plot([p['ra'] for p in points], [p['dec'] for p in points ], 'rx', markersize=2)

    if pointing is not None:
        # print('pointing', pointing['ra'], pointing['dec'])
        ax.plot([pointing['ra']], [pointing['dec']], 'k+')
        ax.text(pointing['ra']+0.05, pointing['dec']+0.06, 'P')

    if title is not None:
        ax.set_title(title)

    ax.set_xlim((82.5, 84.7))
    ax.set_ylim((21.0, 23.1))
    ax.set_aspect('equal')
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    plt.show()

def create_pixel_map(region, pixel_side=0.1):
    for k in ['ra', 'dec', 'rad']:
        if k in region:
            continue
        raise Exception('region data missing {} mandatory key.'.format(k))
    if region['rad'] <= 0:
        raise Exception('region radius must be > 0')
    if pixel_side <= 0:
        raise Exception('pixel side must be > 0')

    region_center = { 'ra': float(region['ra']), 'dec': float(region['dec']) }
    region_rad = utils.get_angle(float(region['rad']))
    pixel_side_angle = utils.get_angle(float(pixel_side))

    # add a 10% to get a bit of margin
    n_pixel_on_diam = 1.1* region_rad / pixel_side_angle
    if n_pixel_on_diam <= 1:
        n_pixel_on_diam = 1

    n_pixel_on_axis = float(math.ceil(n_pixel_on_diam))
    multipliers = np.arange(-1*n_pixel_on_axis, n_pixel_on_axis+1)
    print('region rad: {}, pixel side angle: {}, extimated number of pixel on diam: {}, number of pixel on axis: {}, pixels: {}'.format(region_rad, pixel_side_angle, n_pixel_on_diam, n_pixel_on_axis, len(multipliers)**2))
    pixels_midpoint = []
    for i in multipliers:
        for j in multipliers:
            pixels_midpoint.append({ 'ra':  region_center['ra']  + i * pixel_side_angle.deg,
                                      'dec': region_center['dec'] + j * pixel_side_angle.deg })
    return pixels_midpoint

def select_pixels_in_region(midpoints, region):
    for k in ['ra', 'dec', 'rad']:
        if k in region:
            continue
        raise Exception('region data missing {} mandatory key.'.format(k))
    if region['rad'] <= 0:
        raise Exception('region radius must be > 0')
    if len(midpoints) < 1:
        raise Exception('need at least 1 point to check')
    region_center = utils.get_skycoord(region)
    region_radius = utils.get_angle(region['rad'])
    midpoints_ra = []
    midpoints_dec = []
    for p in midpoints:
        midpoints_ra.append(p['ra'])
        midpoints_dec.append(p['dec'])
    midpoints_coords = SkyCoord(midpoints_ra, midpoints_dec, unit='deg', frame='icrs')
    distances = region_center.separation(midpoints_coords)
    return np.extract(distances < region_radius, midpoints)

def get_thetas(pointing, midpoints):
    for k in ['ra', 'dec']:
        if k in pointing:
            continue
        raise Exception('pointing coord {} is missing.'.format(k))
    if len(midpoints) < 1:
        raise Exception('need at least 1 point to check')
    pnt = utils.get_skycoord(pointing)
    midpoints_ra = []
    midpoints_dec = []
    for p in midpoints:
        midpoints_ra.append(p['ra'])
        midpoints_dec.append(p['dec'])
    midpoints_coords = SkyCoord(midpoints_ra, midpoints_dec, unit='deg', frame='icrs')
    return [ ang.degree for ang in pnt.separation(midpoints_coords) ]
    
def eval_aeff(irf_filename, thetas, energy):
    aeff = EffectiveArea(irf_filename=irf_filename)
    n_points = len(thetas)
    val = 0
    for t in thetas:
        val += aeff.get_aeff_2d_log(t, energy) / n_points
    return val

if __name__ == '__main__':
    irf_filename = sys.argv[1]

    pointing      = { 'ra': 83.6331, 'dec': 22.5145 }
    source_region = { 'ra': 83.6331, 'dec': 22.0145, 'rad': 0.5 }

    for psize in [2.0, 1.0, 0.1, 0.4, 0.5, 0.05, 0.02]:
        # create the pixel map
        pnts = create_pixel_map(source_region, psize)
        plot_region(source_region, pnts, pointing, title='pixels center with {} deg size'.format(psize))

        # select the pixels inside the region
        inside_pixels = select_pixels_in_region(pnts, source_region)
        print('selected pixels: {}'.format(len(inside_pixels)))
        plot_region(source_region, inside_pixels, pointing, title='pixels center with {} deg size'.format(psize))

        # calculate the offsets
        offsets = get_thetas(pointing, inside_pixels)

        # extimate the effective area
        for en in [0.025, 1.0, 150.0]:
            aeff_val = eval_aeff(irf_filename, offsets, en)
            print('Effective area ({0:7.3f} TeV): {1:.5e} cm²'.format(en, aeff_val*1e4))




