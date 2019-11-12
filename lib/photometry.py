from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
import astropy.units as u
import numpy as np
import logging
logging.basicConfig(level=logging.WARN)

# Bintable columns:
# 0  name = 'EVENT_ID'; format = '1J'; bscale = 1; bzero = 2147483648
#    name = 'TIME'; format = '1D'; unit = 's'
#    name = 'RA'; format = '1E'; unit = 'deg'
#    name = 'DEC'; format = '1E'; unit = 'deg'
#    name = 'ENERGY'; format = '1E'; unit = 'TeV'
# 5  name = 'DETX'; format = '1E'; unit = 'deg'
# 6  name = 'DETY'; format = '1E'; unit = 'deg'
#    name = 'MC_ID'; format = '1J'

class Photometrics():
    def __init__(self, args):
        self.events_data = None
        self.events_filename = None
        if 'events_filename' in args:
            self.events_filename = args['events_filename']
            self.events_data = load_data_from_fits_file(self.events_filename)
        elif 'events_list' in args:
            self.events_data = args['events_list']
        self.events_list_checks()

    def events_list_checks(self):
        if self.events_data is None:
            raise Exception('Events data is empy. Need a events list.')
        cols_names = self.events_data.columns.names
        for f in ['RA', 'DEC', 'ENERGY']:
            if f not in cols_names:
                raise Exception("Events data has no '{}' col".format(f))

    @staticmethod
    def load_data_from_fits_file(filename):
        """Load events extension data from a fits file.

        Parameters
        ----------
        filename: str

        Returns
        -------
        events_bintable data
        """
        hdul = fits.open(filename)
        # hdul.info()
        events_bintable = hdul['EVENTS']
        data = events_bintable.data
        hdul.close()
        return data

    def get_skycoord(self, pnt_coord):
        coord = None
        if isinstance(pnt_coord, SkyCoord):
            coord = pnt_coord
        elif isinstance(pnt_coord, dict) and 'ra' in pnt_coord and 'dec' in pnt_coord:
            coord = SkyCoord(ra=pnt_coord['ra'], dec=pnt_coord['dec'], unit='deg', frame='icrs')
        else:
            raise Exception('The input parameter must be a SkyCoord or a { "ra": 12.3, "dec": 45.6 } dictionary.')
        return coord

    def get_angle(self, input_angle):
        ang = None
        if isinstance(input_angle, Angle):
            ang = input_angle
        elif isinstance(input_angle, float):
            ang = Angle(input_angle, unit='deg')
        else:
            raise Exception('The input parameter must be an Angle or a float for decimal degree.')
        return ang

    def region_counter(self, input_center, input_radius, emin=None, emax=None):
        """Counts photons in an input area"""
        region_center = self.get_skycoord(input_center)
        region_radius = self.get_angle(input_radius)

        # filtering...
        condlist = np.full(len(self.events_data.field('ENERGY')), True)
        # ... w/ energy boundaries
        if emin is not None:
            condlist &= self.events_data.field('ENERGY') >= emin
        if emax is not None:
            condlist &= self.events_data.field('ENERGY') <= emax

        events_list = np.extract(condlist, self.events_data)
        # events coordinates from the selected events list
        events_coords = SkyCoord(events_list.field('RA'), events_list.field('DEC'), unit='deg', frame='icrs')
        distances = region_center.separation(events_coords)
        return np.count_nonzero(distances < region_radius)

    def reflected_regions(self, input_pointing_center, input_region_center, input_region_radius):
        """Find regions with reflected algorithm.

        Parameters
        ----------
        input_pointing_center: SkyCoord or dict
        input_region_center: SkyCoord or dict
        input_region_radius: Angle or float

        Returns
        -------
        array of regions
        """
        pointing_center = self.get_skycoord(input_pointing_center)
        region_center = self.get_skycoord(input_region_center)
        region_radius = self.get_angle(input_region_radius)

        # Angular separation of reflected regions. 1.05 factor is to have a margin
        region_diameter = 1.05 * 2.0 * region_radius
        radius = pointing_center.separation(region_center)
        numbers_of_reflected_regions = int(2 * np.pi * radius / region_diameter) # floor down
        regions_offset_angle = Angle(360, unit='deg') / numbers_of_reflected_regions

        regions = []
        # starting from the source region 0, we skip region 1 and region N, so 2..N-1
        starting_pos_angle = pointing_center.position_angle(region_center)
        for i in range(2, numbers_of_reflected_regions-1):
            theta = starting_pos_angle + i * regions_offset_angle
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

    def wobble_regions(self, input_pointing_center, input_region_center, input_region_radius):
        """Return the three background regions starting from pointing and source one.

        Parameters
        ----------
        input_pointing_center: SkyCoord or dict
        input_region_center: SkyCoord or dict
        input_region_radius: Angle or float

        Returns
        -------
        array of regions
        """
        pointing_center = self.get_skycoord(input_pointing_center)
        region_center = self.get_skycoord(input_region_center)
        region_radius = self.get_angle(input_region_radius)
        radius = pointing_center.separation(region_center)
        starting_pos_angle = pointing_center.position_angle(region_center)
        regions = []
        for i in range(1,4):
            theta = starting_pos_angle + i * Angle(90, unit='deg')
            coord_pos = pointing_center.directional_offset_by(theta, radius)
            regions.append({ 'ra': coord_pos.ra.deg, 'dec': coord_pos.dec.deg, 'rad': region_radius.deg })
        return regions

