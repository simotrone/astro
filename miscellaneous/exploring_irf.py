from lib.irf_aeff import IRF
from lib.irf_aeff import EffectiveArea
import sys

irf_filename = sys.argv[1]

irf = IRF(irf_filename)
irf.hdul.info()
eff_area = irf.get_eff_area()
# print('EffArea columns:\n', eff_area.columns)
# print('Effective area:\n', eff_area.data)

bkg = irf.get_extension('BACKGROUND')
print('Background columns:\n', bkg.columns)

# test 1
aeff = EffectiveArea(irf_filename=irf_filename)
# test 2
aeff = EffectiveArea(eff_area_bintable=irf.get_eff_area())
# print(aeff)
print(aeff.columns())
# print(aeff.energies)
print('1D. offset 0 deg,    25 GeV; aeff: {0:.2e} m²'.format(aeff.get_aeff_1d_log(0.0, 0.025)))
print('1D. offset 0.8 deg, 100 TeV; aeff: {0:.2e} m²'.format(aeff.get_aeff_1d_log(0.8, 100.0)))

print('2D. offset 0 deg,    25 GeV; aeff: {0:.2e} m²'.format(aeff.get_aeff_2d_log(0.0, 0.025)))
print('2D. offset 0.8 deg, 100 TeV; aeff: {0:.2e} m²'.format(aeff.get_aeff_2d_log(0.8, 100.0)))
