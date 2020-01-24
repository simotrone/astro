import gammalib
import os
import sys

# all in MeV
def flux_eval_following_ctools(energies, prefactor=5.7e-16, index=-2.48, epivot=0.3e6):
    print("energies", energies)
    _exp = index + 1
    print("_exp", _exp)
    print("pre", prefactor * epivot / _exp)
    print("e1", (energies[1]/epivot)**_exp)
    print("e2", (energies[0]/epivot)**_exp)
    res = prefactor * epivot / _exp * ( (energies[1]/epivot)**_exp - (energies[0]/epivot)**_exp )
    return res

def eflux_eval_following_ctools(energies, prefactor=5.7e-16, index=-2.4, epivot=0.3e6):
    _exp = index + 2
    _epivot2 = epivot**2
    res = prefactor * _epivot2 / _exp * ( (energies[1]/epivot)**_exp - (energies[0]/epivot)**_exp )
    return res

def main():
    mev2erg = 1.6021765e-6 # MeV => erg

    input_model = sys.argv[1] 
    models = gammalib.GModels(input_model)
    print(models)
    first_model = models[0]
    m_spect = first_model.spectral()
    print("TS:", first_model.ts())

    emin = gammalib.GEnergy(25, 'GeV')
    emax = gammalib.GEnergy(150, 'TeV')
    r = m_spect.flux(emin, emax)
    er = m_spect.eflux(emin, emax)
    print('gl Flux: {:.4e} ph/cm²/s'.format(r))
    print('fl EFlux: {:.4e} erg/cm²/s'.format(er))
    for par in ['Prefactor', 'PivotEnergy', 'Index']:
        print(par, m_spect[par].value(), m_spect[par].error())

    # pref = m_spect['Prefactor'].value()
    # print('pref', pref)

    print('my flux: {0:.4e} ph/cm²/s'.format(flux_eval_following_ctools([25*1e3, 150*1e6])))
    my_eflux = eflux_eval_following_ctools([25*1e3, 150*1e6]) #, m_spect['Prefactor'].value())
    print('my eflux: {0:.4e} MeV/cm²/s => {1:.4e} erg/cm²/s'.format(my_eflux, my_eflux*mev2erg))

if __name__ == '__main__':
    main()
