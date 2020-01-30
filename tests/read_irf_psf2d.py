from astropy.io import fits
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.ticker as ticker
import numpy as np
import logging
from lib.irf import PSF

if __name__ == '__main__':
    filename = sys.argv[1]
    psf = PSF(irf_filename=filename)

    for i_vals in [ { 'theta': 0.5, 'energy': 1.0 },
                    { 'theta': 1.5, 'energy': 1.0 },
                    { 'theta': 2.5, 'energy': 1.0 }, ]:
        v = psf.get_psf_values(i_vals['theta'], i_vals['energy'])
        #v_1d_log = psf.get_psf_1d_log(i_vals['theta'], i_vals['energy'])
        print('values =>', v) #,"\n  ", v_1d_log)
        print('delta max:', psf.get_psf_delta_max(i_vals['theta'], i_vals['energy']), '@ theta:', i_vals['theta'], ', en:', i_vals['energy'])

    src = { 'ra': 83.6331, 'dec': 22.0145, 'rad': 0.2 }
    pnt = { 'ra': 83.6331, 'dec': 22.5145 }
    energies = [0.025, 0.050, 0.100, 0.200, 1.0, 100, 150]
    # energies = [0.025, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 1.0, 100, 150]

    for energy in energies:
        flux_rate = psf.eval_region_flux_rate(src, pnt, energy)
        print('flux rate at {:7.3f} TeV: {:.3f}'.format(energy, flux_rate[0]))

    # testing the engine
    engines = {}
    for energy in energies:
        psf_eng = psf.get_psf_engine(src, pnt, energy)
        engines[energy] = psf_eng

    for k, eng in engines.items():
        print('Energy: {:7.3f} psf%: {:.3f}'.format(k, eng(0, np.deg2rad(0.2))[0]))

    # test for radius parts
    limits = np.linspace(0, np.deg2rad(0.2), 10)
    # print(limits)
    for k, eng in engines.items():
        vals = [ eng(limits[i], limits[i+1])[0] for i,l in enumerate(limits[:-1]) ]
        print('Energy: {:7.3f} psf%: {:.3f}'.format(k, np.sum(vals)), vals)
    exit(1)
