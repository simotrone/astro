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
import pprint
import math

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
                    'li_ma_2': [],
                }
            }

        # FIXME excluding all the excess_count < 0
        if float(i['excess_count']) < 0:
            continue

        o[key]['data']['seed'].append(int(i['seed']))
        o[key]['data']['ts'].append(float(i['ts'])) # if float(i['ts']) > 0 else 0)
        o[key]['data']['flux'].append(float(i['flux']))
        o[key]['data']['eflux'].append(float(i['eflux']))
        o[key]['data']['N_on'].append(float(i['on_count']))
        o[key]['data']['N_off'].append(float(i['off_count']))
        o[key]['data']['N_exc'].append(float(i['excess_count']))
        o[key]['data']['alpha'].append(float(i['alpha']))
        o[key]['data']['li_ma'].append(float(i['li_ma']) if i['li_ma'] != '' else 0)

        o[key]['data']['li_ma_2'].append(float(i['li_ma'])**2 if i['li_ma'] != '' else 0)

        if float(i["ts"]) < 0:
            logging.warning("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Negative ts {5:.2f}".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"]), float(i["ts"])))
        elif i["li_ma"] is None:
            logging.warning("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Cannot calculate Li&Ma".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"])))
    return o

def data_count(data, key):
    ret = {}
    for data_name, d in data.items():
        ret[data_name] = {
            'len': len(d['data'][key]),
            'count_0+': 0,
            'count_0-': 0,
            'count_0': 0,
        }
        for i in d['data'][key]:
            if i > 0:
                ret[data_name]['count_0+'] += 1
            elif i < 0:
                ret[data_name]['count_0-'] += 1
            else:
                ret[data_name]['count_0'] += 1
    for k, v in ret.items():
        ret[k]['count_0+%'] = v['count_0+'] / v['len']
        ret[k]['count_0-%'] = v['count_0-'] / v['len']
        ret[k]['count_0%'] = v['count_0'] / v['len']
    return ret

# WARNING: this function augment the input data struct
def data_augmentation(data, bins_number=50):
    fields = [
        { 'name': 'ts',    'dyn_bins': False },
        { 'name': 'N_on',  'dyn_bins': False },
        { 'name': 'N_off', 'dyn_bins': False },
        { 'name': 'N_exc', 'dyn_bins': False },
        { 'name': 'li_ma', 'dyn_bins': False },
        { 'name': 'li_ma_2', 'dyn_bins': True },
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

            counts_hist_normalized = counts_hist / np.sum(counts_hist)
            counts_hist_normalized_error = np.sqrt(counts_hist_normalized) / np.sqrt(np.sum(counts_hist))

            data_stats = array_stats(data_arr_ref)
            d['stats'][f_name] = data_stats

            # starting_parameters = [1., data_stats['mean'], data_stats['stdev']] # A, mu, sigma
            # fit_coeff, pvalue_err = fitting_data(gauss, initial_params=starting_parameters, x=bins_centres, y=counts_hist, verbosity=False, name=data_name)

            d['hist'][f_name] = {
                'n_bins': n_bins,
                'counts': counts_hist,
                'counts_error': np.sqrt(counts_hist),
                'bins_edges': bins_edges,
                'bins_centres': bins_centres,
                'bins_width': bins_width,
                # 'fit_coeff':  fit_coeff,
                # 'pvalue_err': pvalue_err,
                'counts_norm': counts_hist_normalized,
                'counts_norm_error': counts_hist_normalized_error,
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

def print_txt_data_summary(data):
    fields = [
      #   h_format,      v_format,               title,             sub_t
        [   '%15s',        '%15s',            'fs ref',      '==========', ],
        [   '%10s',        '%10s',                'RA',              '==', ],
        [   '%10s',        '%10s',               'Dec',             '===', ],
        [    '%6s',         '%6d',              'tmax',            '====', ],
        [    '%6s',         '%6d',             'seeds',           '=====', ],
        [   '%16s', '%9.2f±%6.2f',                'TS',              '==', ],
        [   '%15s', '%8.2f±%6.2f',              'N_on',            '====', ],
        [   '%15s', '%8.2f±%6.2f',             'N_off',           '=====', ],
        [   '%18s', '%9.2e±%6.2e',               'N_s',             '===', ],
        [   '%13s', '%7.3f±%5.3f',             'Li&Ma',           '=====', ],
        # [   '%18s', '%9.2e±%8.2e',   'flux [ph/cm²/s]', '===============', ],
        # [   '%18s', '%9.2e±%8.2e', 'eflux [erg/cm²/s]', '===============', ],
        # [   '%26s', '%10.2f %7.2f %6.2f', 'TS fitting (A, μ, σ)', '=======', ],
        # [   '%23s', '%10.2f %5.2f %5.2f', 'TS pvalue (A, μ, σ)',  '=======', ],
    ]

    header_fmt = ' '.join([r[0] for r in fields]) # headers format
    values_fmt = ' '.join([r[1] for r in fields]) # values format
    print(header_fmt % tuple([r[2] for r in fields])) # titles
    print(header_fmt % tuple([r[3] for r in fields])) # sub_titles separator
    for d in sorted(data.values(), key=lambda i: (-1*i['tmax'], i['ra'], i['dec'])):
        n_seeds = len(d['data']['seed'])
        # alpha_m = array_stats(d['data']['alpha']) # useless
        # sigma_sign = d['stats']['li_ma']['sigma_significance']
        print(values_fmt % (d['name'], d['ra'], d['dec'], d['tmax'], n_seeds,
            d['stats']['ts']['mean'], d['stats']['ts']['stdev'],
      #      flux_m['mean'], flux_m['stdev'],
      #      eflux_m['mean'], eflux_m['stdev'],
            d['stats']['N_on']['mean'], d['stats']['N_on']['stdev'],
            d['stats']['N_off']['mean'], d['stats']['N_off']['stdev'],
            d['stats']['N_exc']['mean'], d['stats']['N_exc']['stdev'],
            d['stats']['li_ma']['mean'], d['stats']['li_ma']['stdev'], ))

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

def dynamic_bin_number(arr, max_val=None, min_val=None):
    return math.ceil(max(arr))
    n = max(arr)-min(arr)
    # print(min(arr), max(arr), math.ceil(max(arr)))
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
def create_hist(ax, data, xlabel=None, color="blue"):
    bins_centres = data['bins_centres']
    bins_width   = data['bins_width']

    ax.set_yscale('log')
    ax.set_xlabel('{0}'.format(xlabel))

    # print('DEBUG')
    # print('Counts:', data['counts'], np.sum(data['counts']))
    # print('Counts err:', data['counts_error'])
    # print('Counts norm:', data['counts_norm'], np.sum(data['counts_norm']))
    # print('Counts norm err:', data['counts_norm_error'])
    # DEBUG: flag to plot things
    normalized = True
    chi2_plot = True
    pvalue_plot = False
    if normalized:
        ax.errorbar(x=bins_centres, y=data['counts_norm'], xerr=bins_width/2, yerr=data['counts_norm_error'], color=color, label='data', fmt=',', alpha=0.8)
        ax.set_ylabel('normalized counts (log)')
        ax.set_xlabel('TS')

    if normalized and chi2_plot:
        k = 1
        mu = 0
        dist = stats.chi2(k, mu)
        ax.plot(bins_centres[1:], dist.pdf(bins_centres[1:]), color="cyan", linestyle="-.", alpha=0.6, label='Χ²')
        ax.plot(bins_centres[1:], dist.pdf(bins_centres[1:])/2, color="orange", linestyle="--", alpha=0.9, label='Χ²/2')

    if False and normalized is False:
        counts_hist = data['counts']
        counts_error = data['counts_error']
        ax.errorbar(x=bins_centres, y=counts_hist, xerr=bins_width/2, yerr=counts_error, alpha=0.5, color=color, label='data', fmt='+')
        ax.set_ylabel('counts (log)')

    if pvalue_plot:
        k = 1
        mu = 0
        dist = stats.chi2(k, mu)
        pvalues = 1-dist.cdf(bins_centres)
        # print("pvalues:", pvalues)
        # print("pvalues err:", np.sqrt(
        pvalues_err = np.sqrt(bins_centres)/np.sum(data['counts'])
        ax.errorbar(x=bins_centres, y=pvalues, xerr=bins_width/2, yerr=pvalues_err, label='pvalue', fmt='k+')
        ax.set_ylabel('pvalue (log)')
        ax.set_xlabel('h')

        tmp_x = np.arange(len(bins_centres))
        ax.plot(tmp_x[1:], 1-dist.cdf(tmp_x[1:]), color="red", linestyle="-.", alpha=0.6, label='P')

    ax.legend(loc='best') # easy way


def plot_data_summary(data, save_img=False):
    rows_num=1
    cols_num=1
    img_format = 'png'
    for d in list(data.values()):
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=cols_num, figsize=(cols_num*4.5, rows_num*4.5))
        fig.suptitle('{} t_{} sec'.format(d['name'], str(d['tmax'])), va='top', ha='center')

        # create_hist(ax[0], d['hist']['li_ma'], color='blue', xlabel='Li & Ma')
        create_hist(ax, d['hist']['li_ma_2'], color='black', xlabel='TS')

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
        # break # DEBUG
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze data from tsv")
    parser.add_argument('dir', help='the dir with tsv files', nargs='+')
    args = parser.parse_args()

    print(args.dir)
    onlyfiles = []
    for mypath in args.dir:
        for f in os.listdir(mypath):
            fn = os.path.join(mypath, f)
            if not os.path.isfile(fn):
                continue
            onlyfiles.append(fn)

    
    original_data = []
    file_count = 0
    for fn in onlyfiles:
        file_data = csvex.load(fn, header=True, sep='\t')
        file_count += 1
        for fd in file_data:
            original_data.append(fd)

    print("File read: {}".format(file_count), file=sys.stderr)
    print("     data: {}".format(len(original_data)), file=sys.stderr)

    ds = data_summary(original_data)
    if not data_summary_is_ok(ds, pointings=1, time_slots=12, different_seeds=100000):
        exit(1)

    # inplace augment
    data_augmentation(ds, bins_number=25)

    if False:
        pp = pprint.PrettyPrinter(indent=4)
        print("Data count (N_exc):")
        pp.pprint(data_count(ds, 'N_exc'))
        # print("Data count (li_ma):")
        # pp.pprint(data_count(ds, 'li_ma'))
        # print("Data count (ts):")
        # pp.pprint(data_count(ds, 'ts'))
    print_txt_data_summary(ds)
    plot_data_summary(ds, save_img=False)

