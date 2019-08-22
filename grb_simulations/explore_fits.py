from astropy.io import fits
import sys
import csv

class TimeSliceExporter:
    def __init__(self, input_filename):
        self.input_filename = input_filename
        self.hdul = fits.open(input_filename)

    def info(self):
        self.hdul.info()

    # times in secs
    # energies in GeV
    # spectra in ph/cm2/s/GeV
    def export(self):
        times    = self.hdul['TIMES'].data
        energies = self.hdul['ENERGIES'].data
        spectra  = self.hdul['SPECTRA'].data
        for i, tsec in enumerate(times):
            print("{0:2d} {1:15f} sec".format(i, tsec[0]))
            time_slice_filename = "spec_{0:02d}.tsv".format(i)
            self.write_slice_tsv(time_slice_filename, energies, spectra[i])

    def write_slice_tsv(self, output_filename, energies, spectra):
        if not len(energies) == len(spectra):
            raise Exception("Need the same number of elements between energies and spectra")
        with open(output_filename, mode='w', newline="\n") as fh:
            writer = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for i, ene in enumerate(energies):
                writer.writerow([ene[0], spectra[i]])

if __name__ == "__main__":
    input_filename = sys.argv[1]
    exporter = TimeSliceExporter(input_filename)
    exporter.info()
    exporter.export()
    exit(0)
