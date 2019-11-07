import argparse
import json
import os
import sys
from lib.exporter.csv import CSVExporter as csvex

def data_summary_is_ok(data):
    # check size
    if len(data) != 5*7:
        print("Data summary length is {} and should be {} (pointings x time_slots)".format(len(data), 5*7, file=sys.stderr))
        return False
    for k in data:
        # check array data
        for sk in data[k]:
            if type(data[k][sk]) != type([]):
                continue
            if len(data[k][sk]) == 5000:
                continue
            print("not enough data for '{}'".format(k), file=sys.stderr)
            print("  key '{}' has {} values and should be {}".format(sk, len(data[k][sk]), 5000), file=sys.stderr)
            return False
    return True

parser = argparse.ArgumentParser(description="Transform json data to tsv")
parser.add_argument("filename", help="the json file")
parser.add_argument("-f", "--force", help="force file overwriting", default=False, action="store_true")
args = parser.parse_args()

dump_fn = args.filename
if not os.path.isfile(dump_fn):
    raise Exception("cannot read {}".format(dump_fn))

print("Load data from: {}".format(dump_fn), file=sys.stderr)
ds = None
with open(dump_fn, "r") as f:
    ds = json.load(f)
if not data_summary_is_ok(ds):
    exit(1)

# dec_0.5_1800 
# dict_keys(['name', 'tmax', 'ra', 'dec', 'event_file', 'model_file', 'onoff_obs', 'seed_array', 'flux_array', 'eflux_array', 'ts_array', 'sqrt_ts_array', 'li_ma_array', 'N_on_array', 'N_off_array', 'N_s_array', 'alpha_array'])

output = {}
for name, d in ds.items():
    print("Name:", name)
    name = d['name']
    # tmax = d['tmax']
    ra   = d['ra']
    dec  = d['dec']
    # event_file = d['event_file']
    # model_file = d['model_file']
    # onoff_obs  = d['onoff_obs']
    seed_arr   = d['seed_array']
    ts_arr     = d['ts_array']
    N_on_arr   = d['N_on_array']
    N_off_arr  = d['N_off_array']
    alpha_arr  = d['alpha_array']
    li_ma_arr  = d['li_ma_array']
    # sqrt_ts_arr = d['sqrt_ts_array']
    N_s_arr    = d['N_s_array']
    flux_arr   = d['flux_array']
    eflux_arr  = d['eflux_array']

    working_dir = name
    for i, seed in enumerate(seed_arr):
        # if seed != 1: # debug check
        #     continue
        filename = os.path.join(working_dir, 'results_{}.tsv'.format(str(seed)))
        if not args.force and os.path.isfile(filename):
            raise Exception("file {} already exists".format(filename))

        if filename not in output:
            output[filename] = []
        z = {
            'name': name,
            'ra': ra,
            'dec': dec,
            'seed': seed,
            'tmax': d['tmax'],
            'ts': ts_arr[i],
            'on_count': N_on_arr[i],
            'off_count': N_off_arr[i],
            'alpha': alpha_arr[i],
            'li_ma': li_ma_arr[i],
            ####
            'flux': flux_arr[i],
            'eflux': eflux_arr[i],
        }
        output[filename].append(z)

    # debug
    # for k,v in d.items():
    #     if k == "name" or k == "event_file":
    #         print(" ", k, ":", v)
    #     else:
    #         print(" ", k)
    # break

print(len(output))

for fn, data in output.items():
    csvex.save(fn, data, headers=list(data[0].keys()), delimiter="\t")

