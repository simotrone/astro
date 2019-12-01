from astropy.io import fits
import numpy as np
from scipy import interpolate
from lib import utils 

class IRF:
    def __init__(self, filename):
        self.filename = filename
        self.hdul = fits.open(self.filename)

    def get_extension(self, name):
        return self.hdul[name]

    def get_eff_area(self):
        return self.get_extension('EFFECTIVE AREA')

class EffectiveArea:
    def __init__(self, irf_filename=None, eff_area_bintable=None):
        self.irf_filename = None
        self.eff_area = None

        # a sort of cache...
        self.energies = None
        self.thetas = None
        self.aeff_matrix = None

        if irf_filename is not None:
            self.irf_filename = irf_filename
            irf = IRF(self.irf_filename)
            self.eff_area = irf.get_eff_area()
        elif eff_area_bintable is not None:
            self.eff_area = eff_area_bintable

        if self.eff_area is None:
            raise Exception('Need an irf or effective area bintable')

        self.get_data_matrices()

    def columns(self):
        return self.eff_area.columns

    def data(self):
        return self.eff_area.data

    def get_data_matrices(self):
        """ returns aeff data matrix, energy[LO,HI] and theta[LO,HI]
        """
        if self.energies is None or self.thetas is None or self.aeff_matrix is None:
            data = self.data()
        if self.aeff_matrix is None:
            self.aeff_matrix = data.field('EFFAREA')[0]
        if self.energies is None:
            self.energies = np.column_stack((data.field('ENERG_LO')[0], data.field('ENERG_HI')[0]))
        if self.thetas is None:
            self.thetas = np.column_stack((data.field('THETA_LO')[0], data.field('THETA_HI')[0]))
        return self.aeff_matrix, self.energies, self.thetas

    def get_aeff_1d_log(self, offset, energy):
        """
        return effective area in [m²]

        Parameters
          offset: Angle
          energy: [TeV]

        This method does 1D interpolation on energy range, managed as log10.
        Theta offset is not interpolated.
        """
        offset_angle = utils.get_angle(offset)
        aeff_matrix, energy_bins, theta_bins = self.get_data_matrices()

        theta_index = None
        for i, tb in enumerate(theta_bins):
            if offset_angle.degree >= tb[0] and offset_angle.degree < tb[1]:
                theta_index = i
                break
        if theta_index is None:
            raise Exception('Theta offset is out of range ({})'.format(offset))

        # energy interpolation
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        energy_fn  = interpolate.interp1d(x = energy_mid, y = aeff_matrix[theta_index])
        return energy_fn(np.log10(energy))

    def get_aeff_2d_log(self, input_offset, input_energy):
        """
        return effective area in [m²]

        Parameters
          offset: Angle
          energy: [TeV]

        This method does 2D interpolation and manage energy as log10.
        """
        offset_angle = utils.get_angle(input_offset)
        aeff_matrix, energy_bins, theta_bins = self.get_data_matrices()

        # energy interpolation
        theta_mid  = [ (t[0]+t[1])/2 for t in theta_bins ]
        energy_mid = [ (e[0]+e[1])/2 for e in np.log10(energy_bins) ]
        energy_fn  = interpolate.interp2d(x = energy_mid, y = theta_mid, z = aeff_matrix)
        return energy_fn(np.log10(input_energy), input_offset)[0]

