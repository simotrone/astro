import argparse
import collections
import logging
import os
import matplotlib.pyplot as plt
from lib.ctoolswrapper import CToolsWrapper

logging.basicConfig(format='%(asctime)s %(levelname)s:\n%(message)s', level=logging.WARNING)

# PYTHONPATH=.. python test_events_generation.py --model ../crab_simulations/crab.xml --save --dir pippo --tmax 100

parser = argparse.ArgumentParser(description="Create spectrum and plot it")
parser.add_argument("input_file", help="The events fits file")
parser.add_argument("--name", help="source name (descriptive)", default="test")
parser.add_argument("-m", "--model", help="source model for simulation", required=True)
# parser.add_argument("-t", "--tmax", help="the final observations time in seconds", type=int, default=1800)
# parser.add_argument("--ra",  help="simulations right ascension degrees", type=float, default=0)
# parser.add_argument("--dec", help="simulations declination degrees", type=float, default=0)
parser.add_argument("--dir",   help="the savings directory (default: data/)", default="data")
# parser.add_argument("--seed",  help="simulations seed (default: 1)", default=1, type=int)
parser.add_argument("--force", help="force the overwriting", default=False, action="store_true")
parser.add_argument("--save",  help="save the outputs", default=False, action="store_true")
parser.add_argument("-v", "--verbose", action="count", default=0)
args = parser.parse_args()

tw = CToolsWrapper({ 'name': args.name,
                     'ra': 0,
                     'dec': 0,
                     'energy_min': 0.03, 
                     'energy_max': 150.0,
                     #'seed': args.seed,
                     }, verbosity=args.verbose)

working_dir = os.path.join(args.dir)
try:
    os.makedirs(working_dir)
except FileExistsError as e:
    logging.warning("The data dir {} already exists".format(working_dir))

output_filename = os.path.join(working_dir, 'test_spectrum.fits')
log_filename = os.path.join(working_dir, 'test_csspec.log')

spec = tw.csspec_run( input_obs_list = args.input_file,
                      input_models  = args.model,
                      output_file = output_filename,
                      log_file = log_filename,
                      force = args.force,
                      save  = args.save )

fits = spec.spectrum()
# from show_spectrum cscript
print(fits)
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

    print(str(got), ts, flx, e_flx, sep='\t')

# Set upper limit errors
yerr = [0.6 * x for x in ul_flux]

# Plot the spectrum
plt.figure(figsize=(8,5))
plt.loglog()
plt.grid()
plt.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs], fmt='ro')
plt.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs], yerr=yerr, uplims=True, fmt='ro')
plt.xlabel('Energy (TeV)')
plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')
plt.title('{} spectrum'.format(args.name))
plt.show()


