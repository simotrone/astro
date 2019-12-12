from astropy.io import fits
import numpy as np
from scipy import interpolate, integrate
from lib import utils 
from astropy.coordinates import SkyCoord
import math

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

    # this method use an energy range to evaluate the aeff.
    # The energy range is binned and weighted with a powerlaw with index = e_index.
    # Each pixel has a radial weigth and a specific column of energies weight.
    # Lower energies have more weigth than higher.
    # if energy range is small, the effect is trascurable - similar to weighted_value_for_region_single_energy method.
    def weighted_value_for_region(self, region, pointing, input_energies, pixel_size=0.05, e_index=-2.4):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energies: a couple of values in TeV (ex: [ 0.025, 1.0 ])
          pixel_size: a value in degree (default: 0.05)
          e_index: is the powerlaw index (default: -2.4)
        """
        if len(input_energies) != 2:
            raise Exception('need two energies')

        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        log_energies = np.log10(input_energies)
        # N steps for every unit of log energy
        steps = int(np.ceil(log_energies[1]-log_energies[0]) * 10)
        energies = 10**np.linspace(log_energies[0], log_energies[1], steps)
        powerlaw = lambda x: x**e_index
        i_full = integrate.quad(powerlaw, input_energies[0], input_energies[1])
        i_partials = [ integrate.quad(powerlaw, energies[i], energies[i+1]) for i,v in enumerate(energies[:-1]) ]
        i_factor = [ p[0]/i_full[0] for p in i_partials ]
        energies_middle = (energies[1:]+energies[:-1])/2
        n_points = len(offsets)
        val = 0
        for t in offsets:
            for i, en in enumerate(energies_middle):
                val += self.get_aeff_2d_log(t, en) * i_factor[i] / n_points
        return val

    # this method use an energy range to evaluate the aeff. The energy range is
    # binned and every matrix cube (pixel distance * energy bin) have the same
    # weight.
    def weighted_value_for_region_v02(self, region, pointing, input_energies, pixel_size=0.05):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energies: a couple of values in TeV (ex: [ 0.025, 1.0 ])
          pixel_size: a value in degree (default: 0.05)
        """
        if len(input_energies) != 2:
            raise Exception('need two energies')

        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        log_energies = np.log10(input_energies)
        diff = np.ceil(log_energies[1]-log_energies[0])
        steps = int(diff * 10) # N steps for every unit of log energy
        energies = 10**np.linspace(log_energies[0], log_energies[1], steps)
        n_points = len(offsets) * len(energies)
        val = 0
        for t in offsets:
            for en in energies:
                val += self.get_aeff_2d_log(t, en) / n_points
        return val

    # Deprecated 2019-12-11
    # this method is old. just a plain output from an array of points and ONE
    # energy (usually the middle point of the energy range)
    def weighted_value_for_region_single_energy(self, region, pointing, energy, pixel_size=0.05):
        """return effective area value [m²] for a specific region

        Parameters
          region:   { 'ra': ..., 'dec': ..., 'rad': ... }
          pointing: { 'ra': ..., 'dec': ... }
          energy: a value in TeV (ex: 0.025, 1.0, 150.0)
          pixel_size: a value in degree (default: 0.05)
        """
        # create a grid of points
        points = self.create_pixel_map(region, pixel_size)
        # select the points inside the region
        internal_points = self.select_points_in_region(points, region)
        # calculate the offsets
        offsets = self.get_thetas(pointing, internal_points)

        n_points = len(offsets)
        val = 0
        for t in offsets:
            val += self.get_aeff_2d_log(t, energy) / n_points
        return val # m2

    # helpers
    @staticmethod
    def create_pixel_map(region, pixel_side):
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

        # +10% to get a bit of margin
        n_pixel_on_rad = 1.1* region_rad / pixel_side_angle
        if n_pixel_on_rad <= 1:
            n_pixel_on_rad = 1

        n_pixel_on_axis = float(math.ceil(n_pixel_on_rad))
        multipliers = np.arange(-1*n_pixel_on_axis, n_pixel_on_axis+1)
        pixels_midpoint = []
        for i in multipliers:
            for j in multipliers:
                pixels_midpoint.append({ 'ra':  region_center['ra']  + i * pixel_side_angle.deg,
                                         'dec': region_center['dec'] + j * pixel_side_angle.deg })
        return pixels_midpoint

    @staticmethod
    def select_points_in_region(midpoints, region):
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

    @staticmethod
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
    
