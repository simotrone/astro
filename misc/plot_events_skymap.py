import argparse
from astropy.io import fits
import matplotlib.pyplot as plt

def get_events_from_file(input_filename):
    hdul = fits.open(input_filename)
    events_bintable = hdul['EVENTS']
    data = events_bintable.data
    hdul.close()
    return data


def data_selection(events):
    t   = events.field('TIME')
    ra  = events.field('RA')
    dec = events.field('DEC')
    ene = events.field('ENERGY')
    ret = { "time":t, "ra":ra, "dec":dec, "ene":ene }
    return ret


def map_plot(data, ra=None, dec=None, plt_bins=10):
    fig = plt.figure(figsize=(6,6))
    grid = plt.GridSpec(4,4, hspace=0.1, wspace=0.1)
    main_ax  = fig.add_subplot(grid[1:, :3])
    top_ax   = fig.add_subplot(grid[0, :3], sharex=main_ax)
    right_ax = fig.add_subplot(grid[1:, 3], sharey=main_ax)

    main_ax.hist2d(data["ra"], data["dec"], bins=plt_bins) # cmin=10
    # track the coord reference
    if ra:
        main_ax.axvline(ra, color="r", linestyle=":", alpha=0.3)
        top_ax.axvline(ra,  color="r", linestyle=":", alpha=0.4)
    if dec:
        main_ax.axhline(dec,  color="r", linestyle=":", alpha=0.3)
        right_ax.axhline(dec, color="r", linestyle=":", alpha=0.4)

    top_ax.set_title("RA")
    top_ax.hist(data["ra"], bins=plt_bins, alpha=0.5, edgecolor="blue", histtype="stepfilled")
    plt.setp(top_ax.get_xticklabels(), visible=False)

    right_ax.set_title("Dec")
    right_ax.hist(data["dec"], bins=plt_bins, alpha=0.5, edgecolor="blue", histtype="stepfilled", orientation="horizontal")
    plt.setp(right_ax.get_yticklabels(), visible=False)

    plt.show()
    plt.close()


def main(args):
    events = get_events_from_file(args.input_file)
    data = data_selection(events)
    map_plot(data, ra=args.ra, dec=args.dec, plt_bins=30)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="show events skymap")
    parser.add_argument('input_file', help='the input fits file')
    parser.add_argument('-ra', type=float, default=None, help='a Right Ascension reference')
    parser.add_argument('-dec', type=float, default=None, help='a Declination reference')
    args = parser.parse_args()
    main(args)
    exit(0)
