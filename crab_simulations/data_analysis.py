from lib.exporter.csv import CSVExporter as csvex
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats
from scipy.optimize import curve_fit
import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import statistics
import sys

# Example:
# PYTHONPATH=../path/to/lib/ python data_analysis.py *tsv

def data_summary_is_ok(data, pointings=None, time_slots=None, different_seeds=None):
    if len(data) != pointings * time_slots:
        logging.warning("Data summary length is {} and should be {} (pointings x time_slots)".format(len(data), pointings*time_slots))
        return False
    # check array data
    for k in data:
        for sk in data[k]:
            if type(data[k][sk]) != type([]):
                continue
            if len(data[k][sk]) == different_seeds:
                continue
            logging.warning("not enough data for '{}'".format(k))
            logging.warning("  key '{}' has {} values and should be {}".format(sk, len(data[k][sk]), different_seeds))
            return False
    return True

def data_summary(all_data_info):
    o = {}
    for i in all_data_info:
        key = i['name'] +'_'+ str(i['tmax'])
        if key not in o:
            o[key] = {
                'name': i['name'],
                'tmax': int(i['tmax']),
                'ra':   float(i['ra'])  if 'ra'  in i else None,
                'dec':  float(i['dec']) if 'dec' in i else None,
                'data': {
                    'seed':  [],
                    'ts':    [],
                    'index_value': [],
                    'index_error': [],
                    'prefactor_value': [],
                    'prefactor_error': [],
                    'pivot_value': [],
                    'pivot_error': [],
                    'flux':  [],
                    'eflux': [],
                    'N_on':  [],
                    'N_off': [],
                    'N_exc': [],
                    'alpha': [],
                    'li_ma': [],
                }
            }

        o[key]['data']['seed'].append(int(i['seed']))
        o[key]['data']['ts'].append(float(i['ts']))
        o[key]['data']['flux'].append(float(i['flux']))
        o[key]['data']['eflux'].append(float(i['eflux']))
        o[key]['data']['N_on'].append(float(i['on_count']))
        o[key]['data']['N_off'].append(float(i['off_count']))
        o[key]['data']['N_exc'].append(float(i['excess_count']))
        o[key]['data']['alpha'].append(float(i['alpha']))
        o[key]['data']['li_ma'].append(float(i['li_ma']) if i['li_ma'] != '' else 0)

        if float(i["ts"]) < 0:
            logging.warning("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Negative ts {5:.2f}".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"]), float(i["ts"])))
        elif i["li_ma"] is None:
            logging.warning("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Cannot calculate Li&Ma".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"])))
    return o

# WARNING: this function augment the input data struct
def data_augmentation(data, bins_number=50):
    fields = [
        { 'name': 'ts',    'dyn_bins': False },
        { 'name': 'flux',  'dyn_bins': False },
        { 'name': 'eflux', 'dyn_bins': False },
    ]
    for data_name, d in data.items():
        logging.warning(data_name)
        if 'hist' not in d:
            d['hist'] = {}
        if 'stats' not in d:
            d['stats'] = {}

        for f in fields:
            f_name = f['name']
            data_arr_ref = d['data'][f_name]
            n_bins = dynamic_bin_number(data_arr_ref) if f['dyn_bins'] else bins_number

            # counts histogram
            counts_hist, bins_edges, bin_index_not_used = stats.binned_statistic(data_arr_ref, data_arr_ref, statistic='count', bins=n_bins)
            bins_width = np.array(np.diff(bins_edges), float)
            bins_centres = (bins_edges[:-1] + bins_edges[1:])/2

            # counts_hist_normalized = counts_hist / bins_width / np.sum(counts_hist)

            data_stats = array_stats(data_arr_ref)
            d['stats'][f_name] = data_stats

            starting_parameters = [1., data_stats['mean'], data_stats['stdev']] # A, mu, sigma
            fit_coeff, pvalue_err = fitting_data(gauss, initial_params=starting_parameters, x=bins_centres, y=counts_hist, verbosity=False, name=data_name)

            d['hist'][f_name] = {
                'n_bins': n_bins,
                'counts': counts_hist,
                'bins_edges': bins_edges,
                'bins_centres': bins_centres,
                'bins_width': bins_width,
                'fit_coeff':  fit_coeff,
                'pvalue_err': pvalue_err,
                # 'counts_norm': counts_hist_normalized,
                # 'fit_coeff_norm': fit_coeff_norm,
                # 'pvalue_err_norm': pvalue_err_norm,
            }
    return data

def array_stats(arr):
    stat = {
        "n": len(arr),
        "mean": statistics.mean(arr),
        "stdev": statistics.pstdev(arr),
        "median": statistics.median(arr),
    }
    return stat

def print_data_summary(data):
    fields = [
      #   h_format,      v_format,               title,             sub_t
        [   "%15s",        "%15s",            "fs ref",      "==========", ],
        [   "%10s",        "%10s",                "RA",              "==", ],
        [   "%10s",        "%10s",               "Dec",             "===", ],
        [    "%6s",         "%6d",              "tmax",            "====", ],
        [    "%6s",         "%6d",             "seeds",           "=====", ],
        [   "%16s", "%9.2f±%6.2f",                "TS",              "==", ],
        [   "%18s", "%9.2e±%8.2e",   "flux [ph/cm²/s]", "===============", ],
        [   "%18s", "%9.2e±%8.2e", "eflux [erg/cm²/s]", "===============", ],
        [   '%26s', '%10.2f %7.2f %6.2f', 'TS fitting (A, μ, σ)', '=======', ],
        [   '%23s', '%10.2f %5.2f %5.2f', 'TS pvalue (A, μ, σ)',  '=======', ],
    ]

    header_fmt = " ".join([r[0] for r in fields]) # headers format
    values_fmt = " ".join([r[1] for r in fields]) # values format
    print(header_fmt % tuple([r[2] for r in fields])) # titles
    print(header_fmt % tuple([r[3] for r in fields])) # sub_titles separator
    for d in sorted(data, key=lambda i: (-1*i["tmax"], i["ra"], i["dec"])):
        n_seeds = len(d['data']['seed'])
        ts_m    = array_stats(d['data']['ts'])
        flux_m  = array_stats(d['data']['flux'])
        eflux_m = array_stats(d['data']['eflux'])
        # print(d)
        print(values_fmt % (d["name"], d["ra"], d["dec"], d["tmax"], n_seeds,
            ts_m["mean"], ts_m["stdev"],
            flux_m["mean"], flux_m["stdev"],
            eflux_m["mean"], eflux_m["stdev"],
            d['hist']['ts']['fit_coeff'][0],
            d['hist']['ts']['fit_coeff'][1],
            abs(d['hist']['ts']['fit_coeff'][2]),
            d['hist']['ts']['pvalue_err'][0],
            d['hist']['ts']['pvalue_err'][1],
            d['hist']['ts']['pvalue_err'][2] ))

def fitting_data(curve_fn, initial_params=[], x=[], y=[], verbosity=False, name=None):
    res = curve_fit(curve_fn, x, y, p0=initial_params, full_output=verbosity)
    coeff, var_matrix = res[:2]
    if (len(res) > 2):
        infodict, errmsg, ier = res[2:]
        print("infodict: {}\nerrmsg: {}\nier: {}".format(infodict, errmsg, ier))
    perr = np.sqrt(np.diag(var_matrix))
    print("Curve fit params: {}".format(name))
    print("{0:>10s}  {1:9s}  {2:9s}".format("param no.", "value", "error"))
    for i, c in enumerate(coeff):
        print("{0:10d}  {1:+8.6e}  {2:+8.6e}".format(i, c, perr[i]))
    return coeff, perr

def gauss(x, *params):
    A, mu, sigma = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * np.exp(exp_num / exp_den)

def gauss0(x, *params):
    A, mu, sigma = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * 1. / (2. * np.pi * sigma**2)* np.exp(exp_num / exp_den)

# no good with the quote
def gauss2(x, *params):
    A, mu, sigma, Q = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * 1. / (2. * np.pi * sigma**2)* np.exp(exp_num / exp_den) + Q

def dynamic_bin_number(arr, max_val=None, min_val=None):
    n = max(arr)-min(arr)
    if max_val is not None and n > max_val:
        n = max_val
    if min_val is not None and n < min_val:
        n = min_val
    if n < 1:
        n = 2
    return int(n)

# seaborn graph with distplot. Same data, same gaussian loc/scale
#   import seaborn as sns, numpy as np
#   print(array_stats(d["ts_array"]))
#   print(norm.fit(d["ts_array"]))
#   sns.distplot(d["ts_array"], bins=50, kde=True, fit=norm, norm_hist=False)# , norm_hist=True) #array, bins=n_bins, fit=norm, norm_hist=True
def create_hist(ax, data, data_stats, xlabel=None, color="blue"):
    bins_centres = data['bins_centres']
    bins_width   = data['bins_width']
    fit_params  = data['fit_coeff']
    fitted_hist = gauss(bins_centres, *fit_params)
    counts_hist = data['counts']

    ax.bar(bins_centres, height=counts_hist, width=bins_width, alpha=0.5, edgecolor=color, color=color, label='data')

    # ax.plot(bins_centres, stats.norm.pdf(bins_centres, data_stats["mean"], data_stats["stdev"]), color="orange", linestyle="--", alpha=0.9, label='stats.norm\nμ:{0:.2e}\nσ:{1:.2e}'.format(data_stats['mean'], data_stats['stdev']))

    # gauss fit
    ax.plot(bins_centres, fitted_hist, linestyle="-.", color="green")

    ax.legend( [ Patch(facecolor=color, edgecolor=color),
                 Line2D([0],[0], color="green", linestyle="-.") ],
               [ 'Data. bins:{0}\nμ:{1:.2e}\nσ:{2:.2e}'.format(data['n_bins'], data_stats['mean'], data_stats['stdev']),
                 'Curve fit\nμ:{0:.2e}\nσ:{1:.2e}'.format(fit_params[1], abs(fit_params[2])) ],
               loc='best' )

    ax.set_xlabel('{0}'.format(xlabel))
    ax.set_ylabel('counts')

def plot_data_summary(data, save_img=False):
    rows_num=1
    cols_num=3
    img_format = 'png'
    for d in data:
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=cols_num, figsize=(cols_num*4.5, rows_num*4.5))
        fig.suptitle('{} t_{} sec'.format(d['name'], str(d['tmax'])), va='top', ha='center')

        create_hist(ax[0], d['hist']['ts'],    d['stats']['ts'],    color='blue',  xlabel='Model TS')
        create_hist(ax[1], d['hist']['flux'],  d['stats']['flux'],  color='cyan',  xlabel='flux [ph/cm²/s]')
        create_hist(ax[2], d['hist']['eflux'], d['stats']['eflux'], color='green', xlabel='eflux [erg/cm²/s]')

        # Important: first tight_layout(), after adjust for the title
        fig.tight_layout()
        fig.subplots_adjust(top=0.90)
        if save_img:
            img_filename = "{0}_{1:04d}.{2}".format(d["name"], d["tmax"], img_format)
            plt.savefig(img_filename, format=img_format)
            plt.close()
            logging.debug("saving {}".format(img_filename))
        else:
            plt.show()
            plt.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze data from tsv")
    parser.add_argument('files', help='the tsv file', nargs='+')
    args = parser.parse_args()
    
    original_data = []
    file_count = 0
    for fn in args.files:
        file_data = csvex.load(fn, header=True, sep='\t')
        file_count += 1
        for fd in file_data:
            original_data.append(fd)

    print("File read: {}".format(file_count), file=sys.stderr)
    print("     data: {}".format(len(original_data)), file=sys.stderr)

    ds = data_summary(original_data)
    if not data_summary_is_ok(ds, pointings=5, time_slots=8, different_seeds=5000):
        exit(1)

    # inplace augment
    data_augmentation(ds, bins_number=50)
    print_data_summary(ds.values())
    plot_data_summary(list(ds.values()), save_img=True)

