import argparse
from exporter.timeslice import TimeSliceExporter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an obeservations list from a spectra fits")
    parser.add_argument("input_fits", help="the fits file with times/energies/spectra data")
    parser.add_argument("-m", "--model", help="the xml template")
    parser.add_argument("-d", "--dir", help="the savings directory (default: data/)", default="data")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("--tmax", help="the final observations time in seconds", type=int)
    args = parser.parse_args()

    exporter = TimeSliceExporter(args.input_fits, model_template=args.model, savings_dir=args.dir, verbosity=args.verbose, tmax=args.tmax)
    slices = exporter.export()
    print(slices)
    exit(0)
