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

    def get_region_center(self, input_center):
        region_center = None
        if isinstance(input_center, SkyCoord):
            region_center = input_center
        elif isinstance(input_center, dict) and 'ra' in input_center and 'dec' in input_center:
            region_center = SkyCoord(ra=input_center['ra'], dec=input_center['dec'], unit='deg', frame='icrs')
        else:
            raise Exception('The region center must be a SkyCoord or a { "ra": 12.3, "dec": 45.6 } dictionary.')
        return region_center

    def get_region_radius(self, input_radius):
        region_radius = None
        if isinstance(input_radius, Angle):
            region_radius = input_radius
        elif isinstance(input_radius, float):
            region_radius = Angle(input_radius, unit='deg')
        else:
            raise Exception('The region radius must be an Angle or a float for decimal degree.')
        return region_radius

    def region_counter(self, input_center, input_radius, emin=None, emax=None):
        region_center = self.get_region_center(input_center)
        region_radius = self.get_region_radius(input_radius)

        # filtering...
        condlist = [True] * len(self.events_data.field('ENERGY'))
        # ... w/ energy boundaries
        if emin is not None:
            condlist &= self.events_data.field('ENERGY') >= emin
        if emax is not None:
            condlist &= self.events_data.field('ENERGY') <= emax

        events_list = np.extract(condlist, self.events_data)
        coord = SkyCoord(events_list.field('RA'), events_list.field('DEC'), unit='deg', frame='icrs')

        distances = region_center.separation(coord)
        return np.count_nonzero(distances < region_radius)

