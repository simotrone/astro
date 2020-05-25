import argparse
from astropy.io import fits

def show_header_cards(hdu):
    print('HDU cards:')
    for c in hdu.header.cards:
        print('   ',c)
    
def show_header(hdu):
    print('HDU extension: {}'.format(hdu.header.get('XTENSION', '-')))
    if args.cards:
        show_header_cards(hdu)
        return

    if hasattr(hdu, 'columns'):
        cols = hdu.columns
        fmt = '{:>10s}: '
        for i in cols.names:
            fmt += ' {:>10s}'
        print('HDU Table:')
        print(fmt.format('name', *cols.names))
        print(fmt.format('format', *cols.formats))
        print(fmt.format('unit', *cols.units))
    else:
        print('HDU Table: -')

def show_data(hdu):
    """
    hdu.data is a FITS_rec
    https://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.FITS_rec
    """
    if hdu.data is None:
        print('HDU Data: -')
        return

    # print(hdu.data)
    # hdu.data.names == hdu.columns.names
    aggr = { 'min': [], 'max': [], 'mean': [] }
    fmt_s = '{:>10s}: '
    fmt_v = '{:>10s}: '
    for n in hdu.data.names:
        vals = hdu.data[n]
        aggr['min'].append(vals.min())
        aggr['max'].append(vals.max())
        aggr['mean'].append(vals.mean())
        fmt_s += ' {:>10s}'
        fmt_v += ' {:>10.2e}'
    print('HDU Data:')
    print(fmt_s.format('name', *hdu.data.names))
    print(fmt_v.format('min', *aggr['min']))
    print(fmt_v.format('max', *aggr['max']))
    print(fmt_v.format('mean', *aggr['mean']))
        

def show_hdu(hdu):
    "Show Header Data Unit (HDU)"
    print(f'HDU name: {hdu.name}')
    show_header(hdu)
    show_data(hdu)

def main(args):
    hdul = fits.open(args.input_file)
    hdul.info()

    if args.header:
        show_hdu(hdul[args.header])

    exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="inspect a fits file")
    parser.add_argument('input_file', help='the input fits file')
    # parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('-H', '--header', default=None, help='header name (case insensitive)')
    parser.add_argument('-c', '--cards', default=None, action='store_true', help='show all cards in header')
    args = parser.parse_args()
    main(args)
