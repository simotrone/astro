from lib.photometry import Photometrics 
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def load_fits(filename):
    hdulist = fits.open(filename)
    hdulist.info()
    data = hdulist[0].data
    return data

def plot_regions(regions, pointing=None, title=None, data=None):
    fig, ax = plt.subplots()
    for r in regions:
        rad   = r['rad']   if 'rad'   in r else 0.2
        color = r['color'] if 'color' in r else 'r'
        circle = plt.Circle((r['ra'], r['dec']), rad, color=color, linestyle='--', fill=False)
        ax.add_artist(circle)
        if color == 'green':
            #ax.text(r['ra']-rad, r['dec']+rad+0.05, 'Source')
            ax.text(r['ra']-0.13, r['dec']+0.05, 'Source', fontsize=10)
            ax.plot([r['ra']], [r['dec']], 'k*')
    if pointing is not None:
        ax.plot([pointing['ra']], [pointing['dec']], 'k+')
        ax.text(pointing['ra']-0.05, pointing['dec']+0.05, 'P')
        if False: # arrow
            ax.arrow(x=pointing['ra'], y=pointing['dec'], dx=0.55, dy=0, head_width=0.05, head_length=0.05, ls='--')
            ax.text(pointing['ra']+0.15, pointing['dec']-0.10, '0.4-0.8°')
    if title is not None:
        ax.set_title(title)
    if data is not None:
        print(data)
        ax.imshow(data, cmap=plt.cm.viridis)
    ax.set_xlim((pointing['ra']-1.0, pointing['ra']+1.0))
    ax.set_ylim((pointing['dec']-1.0, pointing['dec']+1.0))
    ax.set_aspect('equal')
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.invert_xaxis()
    plt.show()

def create_map(source, ra_distance=0.5, dec_distance=0.5):
    deltas = [ { 'ra':    ra_distance, 'dec':         0.0 },
               { 'ra':            0.0, 'dec': -1*dec_distance },
               { 'ra': -1*ra_distance, 'dec':         0.0 },
               { 'ra':            0.0, 'dec':    dec_distance }, ]
    pointings = []
    for i, d in enumerate(deltas):
        pointings.append({ 'ra': source['ra']+d['ra'], 'dec': source['dec']+d['dec'], 'name': 't{}'.format(i) })

    for pnt in pointings:
        wobble_regions = Photometrics.wobble_regions(pnt, source, 0.2)
        print('=> pointing:', pnt, '\n   source:', source, '\n   wobble regions:', wobble_regions)
        wobble_regions.append(source)
        # plot_regions(wobble_regions, pointing=pnt, title=pnt['name'], data=fits_data[0]['data'])
        plot_regions(wobble_regions, pointing=pnt, title='observation {}'.format(pnt['name']))

if __name__ == "__main__":
    files = sys.argv[1:]
    if False: # TODO skymap under the hood 
        print(files)
        fits_data = []
        for f in files:
            fits_data.append({ 'file': f, 'data': np.flip(load_fits(f)) })
        fig, ax = plt.subplots()
        ax.imshow(fits_data[0]['data'], cmap=plt.cm.viridis)
        circle = plt.Circle((20, 20), 1)
        #, linestyle='--', fill=False)
        ax.add_artist(circle)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        # plt.colorbar()
        plt.show()
        exit(1)

    crab_source = { 'ra': 83.63,  'dec': 22.01,   'color': 'green' }
    grb_source  = { 'ra': 33.057, 'dec': -51.841, 'color': 'green' }
    create_map(crab_source, ra_distance=0.5, dec_distance=0.5)
    create_map(grb_source,  ra_distance=0.8, dec_distance=0.5)
