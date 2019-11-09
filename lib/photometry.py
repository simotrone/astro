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
        region_center = self.get_skycoord(input_center)
        region_radius = self.get_angle(input_radius)

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

