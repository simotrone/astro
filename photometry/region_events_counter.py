import argparse
from lib.photometry import Photometrics

def main(args):
    phm = Photometrics({ 'events_filename': args.input_file })
    region = {
        'ra': args.right_ascension,
        'dec': args.declination,
    }
    region_count = phm.region_counter(region, args.region_radius)
    print('File: ', args.input_file)
    print('Region center ', region, 'with radius', args.region_radius, 'deg')
    print('Photon count: ', region_count)

    exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="count events in a specific region from a events fits file")
    parser.add_argument('input_file', help='the input fits file')
    parser.add_argument('-ra', '--right-ascension', help='the source right ascension (deg)', type=float)
    parser.add_argument('-dec', '--declination', help='the source declination (deg)', type=float)
    parser.add_argument('-rad', '--region-radius', help='the region radius (default: 0.2°)', default=0.2, type=float)
    args = parser.parse_args()
    main(args)

