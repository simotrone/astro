import argparse
import collections
import logging
import os
import matplotlib.pyplot as plt
from lib.ctoolswrapper import CToolsWrapper
import gammalib

logging.basicConfig(format='%(asctime)s %(levelname)s:\n%(message)s', level=logging.WARNING)

# PYTHONPATH=.. python test_events_generation.py --model ../crab_simulations/crab.xml --save --dir pippo --tmax 100
def read_spectrum_fits(fits):
    table = fits.table(1)
    c_energy = table['Energy']
    c_ed     = table['ed_Energy']
    c_eu     = table['eu_Energy']
    c_flux   = table['Flux']
    c_eflux  = table['e_Flux']
    c_ts     = table['TS']
    c_upper  = table['UpperLimit']
    
    # Initialise arrays to be filled
    energies    = []
    flux        = []
    ed_engs     = []
    eu_engs     = []
    e_flux      = []
    ul_energies = []
    ul_ed_engs  = []
    ul_eu_engs  = []
    ul_flux     = []
    
    # Loop over rows of the file
    for i in range(table.nrows()):
        # Get Test Statistic, flux and flux error
        ts    = c_ts.real(i)
        flx   = c_flux.real(i)
        e_flx = c_eflux.real(i)
    
        # If Test Statistic is larger than 9 and flux error is smaller than flux then append flux plots ...
        got = None
        if ts > 9.0 and e_flx < flx:
            got = True
            energies.append(c_energy.real(i))
            flux.append(c_flux.real(i))
            ed_engs.append(c_ed.real(i))
            eu_engs.append(c_eu.real(i))
            e_flux.append(c_eflux.real(i))
        # ... otherwise append upper limit
        else:
            got = False
            ul_energies.append(c_energy.real(i))
            ul_flux.append(c_upper.real(i))
            ul_ed_engs.append(c_ed.real(i))
            ul_eu_engs.append(c_eu.real(i))
    
        # logging.warning(str(got), ts, flx, e_flx, sep='\t')
    return { 'energies': energies,
             'flux':     flux,
             'ed_engs':  ed_engs,
             'eu_engs':  eu_engs,
             'e_flux':   e_flux,
             'ul_energies': ul_energies,
             'ul_ed_engs':  ul_ed_engs,
             'ul_eu_engs':  ul_eu_engs,
             'ul_flux':     ul_flux, }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create spectrum from a lot of csspec fits files and plot them")
    parser.add_argument("input_files", help="The csspec fits file", nargs='+')
    # parser.add_argument("--save",  help="save the outputs", default=False, action="store_true")
    args = parser.parse_args()
    
    data = []
    for fn in args.input_files:
        spectrum_fits = gammalib.GFits(fn)
        # logging.warning(spectrum_fits)
        data.append(read_spectrum_fits(spectrum_fits))
    
    # Set upper limit errors
    
    # Plot the spectra
    plt.figure(figsize=(8,5))
    plt.loglog()
    plt.grid()

    plt.errorbar(data[0]['energies'], data[0]['flux'], yerr=data[0]['e_flux'], xerr=[data[0]['ed_engs'], data[0]['eu_engs']], fmt='ro', label="crab")
    plt.errorbar(data[0]['ul_energies'], data[0]['ul_flux'], xerr=[data[0]['ul_ed_engs'], data[0]['ul_eu_engs']], yerr=[0.6 * x for x in data[0]['ul_flux']], uplims=True, fmt='ro')

    plt.errorbar(data[1]['energies'], data[1]['flux'], yerr=data[1]['e_flux'], xerr=[data[1]['ed_engs'], data[1]['eu_engs']], fmt='bo', label="grb afterflow")
    plt.errorbar(data[1]['ul_energies'], data[1]['ul_flux'], xerr=[data[1]['ul_ed_engs'], data[1]['ul_eu_engs']], yerr=[0.6 * x for x in data[1]['ul_flux']], uplims=True, fmt='bo')

    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')
    plt.legend()
    # plt.title('{} spectrum'.format(args.name))
    plt.show()


    exit(0)

# tw = CToolsWrapper({ 'name': args.name,
#                      'ra': 0,
#                      'dec': 0,
#                      'energy_min': 0.03, 
#                      'energy_max': 150.0,
#                      #'seed': args.seed,
#                      }, verbosity=args.verbose)
# 
# working_dir = os.path.join(args.dir)
# try:
#     os.makedirs(working_dir)
# except FileExistsError as e:
#     logging.warning("The data dir {} already exists".format(working_dir))
# 
# output_filename = os.path.join(working_dir, 'test_spectrum.fits')
# log_filename = os.path.join(working_dir, 'test_csspec.log')
# 
# spec = tw.csspec_run( input_obs_list = args.input_file,
#                       input_models  = args.model,
#                       output_file = output_filename,
#                       log_file = log_filename,
#                       force = args.force,
#                       save  = args.save )
# 
# fits = spec.spectrum()
# from show_spectrum cscript

