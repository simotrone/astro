from astropy.io import fits
import sys
import csv

def write_time_slice_tsv(filename, energies, spectra):
    if not len(energies) == len(spectra):
        raise Exception("Need same number of elements among energies and spectra")
    with open(filename, mode='w') as fh:
        writer = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, ene in enumerate(energies):
            writer.writerow([ene[0], spectra[i]])

# times in secs
# energies in GeV
# spectra in ph/cm2/s/GeV
def get_data_from_file(input_fn):
    hdul = fits.open(input_filename)
    hdul.info() # DEBUG
    return (hdul['TIMES'].data, hdul['ENERGIES'].data, hdul['SPECTRA'].data)

if __name__ == "__main__":
    input_filename = sys.argv[1]
    times, energies, spectra = get_data_from_file(input_filename)
    for i, tsec in enumerate(times):
        print("%2d %15f sec" % (i, tsec[0]))
        time_slice_filename = "spec_{0:02d}.tsv".format(i)
        write_time_slice_tsv(time_slice_filename, energies, spectra[i])
    exit(0)
