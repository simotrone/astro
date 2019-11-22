import gammalib as g
import sys
import numpy as np

gobs = g.GCTAObservation()
caldb = g.GCaldb('cta', 'prod3b')
irf = 'South_z20_average_30m'
gobs.response(irf, caldb)

rsp = gobs.response()
# print(rsp)
"""
=== GCTAResponseIrf ===
 Caldb mission .............: cta
 Caldb instrument ..........: prod3b
 Response name .............: South_z20_average_30m
 Energy dispersion .........: Not used
 Safe energy range .........: undefined
=== GCaldb ===
 Database root .............: /home/sim/anaconda3/envs/cta_pipe/share/caldb
 Selected Mission ..........: CTA
 Selected Instrument .......: PROD3B
 Calibration Index File ....: /home/sim/anaconda3/envs/cta_pipe/share/caldb/data/cta/prod3b/caldb.indx
 Number of entries .........: 144
=== GCTAAeff2D ===
 Filename ..................: /home/sim/anaconda3/envs/cta_pipe/share/caldb/data/cta/prod3b-v1/bcf/South_z20_average_30m/irf_file.fits
 Number of energy bins .....: 42
 Number of offset bins .....: 6
 Log10(Energy) range .......: 0.0125892544165254 - 199.526229858398 TeV
 Offset angle range ........: 0 - 6 deg
 Lower energy threshold ....: not specified
 Upper energy threshold ....: not specified
 Radius cut ................: none
=== GCTAPsf2D ===
 Filename ..................: /home/sim/anaconda3/envs/cta_pipe/share/caldb/data/cta/prod3b-v1/bcf/South_z20_average_30m/irf_file.fits
 Number of energy bins .....: 25
 Number of offset bins .....: 6
 Log10(Energy) range .......: 0.0125892544165254 - 1258.92541503906 TeV
 Offset angle range ........: 0 - 6 deg
=== GCTAEdisp2D ===
 Filename ..................: /home/sim/anaconda3/envs/cta_pipe/share/caldb/data/cta/prod3b-v1/bcf/South_z20_average_30m/irf_file.fits
 Number of energy bins .....: 500
 Number of migration bins ..: 300
 Number of offset bins .....: 6
 Log10(Energy) range .......: 0.00501187238842249 - 501.187225341797 TeV
 Migration range ...........: 0 - 3
 Offset angle range ........: 0 - 6 deg
=== GCTABackground3D ===
 Filename ..................: /home/sim/anaconda3/envs/cta_pipe/share/caldb/data/cta/prod3b-v1/bcf/South_z20_average_30m/irf_file.fits
 Number of DETX bins .......: 36
 Number of DETY bins .......: 36
 Number of energy bins .....: 21
 DETX range ................: -6 - 6 deg
 DETX range ................: -6 - 6 deg
 Energy range ..............: 0.0125892544165254 - 199.526229858398 TeV
=== GCTAResponseCache ===
 Number of cached values ...: 0
=== GCTAResponseCache ===
 Number of cached values ...: 0
 """

aeff = rsp.aeff()
print(aeff)
print(aeff.table())
nbins = 10
emin = g.GEnergy(0.03, 'TeV')
emax = g.GEnergy(150.0, 'TeV')

def compute_effective_area_for_energy_interval(emin, emax, nbins=10):
    log_emin = emin.log10TeV()
    delta_e = emax.log10TeV() - log_emin
    log_e_bin = delta_e/nbins
    print('ene range: {}-{}\nEmin log10Tev: {}; delta: {}; Ebin: {}'.format(emin, emax, log_emin, delta_e, log_e_bin))
    ref_min = g.GEnergy()
    ref_max = g.GEnergy()
    for n in np.arange(nbins):
        ref_min.log10(log_emin, 'TeV')
        ref_max.log10(log_emin+log_e_bin, 'TeV')
        area = compute_effective_area_at_minimum_and_max_energy(ref_min, ref_max, 10)
        print('{0}  => {1:6.2f}-{2:6.2f} => {3:.2e} cm²'.format(n, ref_min.TeV(), ref_max.TeV(), area))
        log_emin = ref_max.log10TeV()

# see code in ctobssim#1314
def compute_effective_area_at_minimum_and_max_energy(emin, emax, nbins=10):
    logE = emin.log10TeV()
    deltaE = emax.log10TeV() - logE
    logEbin = deltaE/nbins
    # print('energy range: {}-{}\nEmin log10Tev: {}; delta: {}; Ebin: {}'.format(emin, emax, logE, deltaE, logEbin))
    area = 0
    ref = g.GEnergy()
    for n in np.arange(nbins+1):
        a_max = aeff.max(logE, 0.0, 0.0)
        ref.log10TeV(logE) # ref.log10(logE, 'TeV')
        # print(" {0:2d} en: {1:5.2f} {2:6.2f} TeV, Aeff max: {3:.2e}".format(n, logE, ref.TeV(), a_max))
        logE += logEbin
        if a_max > area:
            area = a_max
    return area*2

area = compute_effective_area_at_minimum_and_max_energy(emin, emax, nbins)
print("Final max area: {0:.2e} m² for {1}-{2}".format(area, emin, emax))
compute_effective_area_for_energy_interval(emin, emax, nbins)
# irf = rsp.irf()
# print(irf)
