import argparse
import os
import sys
import statistics
from lib.exporter.csv import CSVExporter as csvex
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np

# Example:
# PYTHONPATH=../path/to/lib/ python data_analysis.py *tsv

def data_summary_is_ok(data):
    if len(data) != 5*7:
        print("Data summary length is {} and should be {} (pointings x time_slots)".format(len(data), 5*7), file=sys.stderr)
        return False
    for k in data:
        # check array data
        for sk in data[k]:
            if type(data[k][sk]) != type([]):
                continue
            if len(data[k][sk]) == 5000:
                continue
            print("not enough data for '{}'".format(k), file=sys.stderr)
            print("  key '{}' has {} values and should be {}".format(sk, len(data[k][sk]), 5000))
            return False
    return True

def data_summary(all_data_info):
    o = {}
    for i in all_data_info:
        key = i["name"] +"_"+ str(i["tmax"])
        if key not in o:
            o[key] = {
                "name": i["name"],
                "tmax": int(i["tmax"]),
                "ra":   float(i["ra"]),
                "dec":  float(i["dec"]),
                "seed_array":  [],
                "flux_array":  [],
                "eflux_array": [],
                "ts_array":    [],
                "li_ma_array": [],
                "N_on_array":  [],
                "N_off_array": [],
                "N_s_array":   [],
                "alpha_array": [],
            }

        o[key]["seed_array"].append(int(i["seed"]))
        o[key]["flux_array"].append(float(i["flux"]))
        o[key]["eflux_array"].append(float(i["eflux"]))
        o[key]["ts_array"].append(float(i["ts"]))
        o[key]["li_ma_array"].append(float(i["li_ma"]) if i["li_ma"] is not None else 0)
        o[key]["N_on_array"].append(float(i["on_count"]))
        o[key]["N_off_array"].append(float(i["off_count"]))
        o[key]["alpha_array"].append(float(i["alpha"]))

        N_s = float(i['on_count']) - float(i['alpha']) * float(i['off_count'])
        o[key]["N_s_array"].append(N_s)

        if float(i["ts"]) < 0:
            print("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Negative ts {5:.2f}".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"]), float(i["ts"])), file=sys.stderr)
        elif i["li_ma"] is None:
            print("{0:15s} ({1:.0f} on, {2:2.0f} off, {3:3d} seed, {4:4d} tmax): Cannot calculate Li&Ma".format(i["name"], float(i["on_count"]), float(i["off_count"]), int(i["seed"]), int(i["tmax"])), file=sys.stderr)
    return o

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
      #   h_format,      v_format,             title,             sub_t
        [   "%15s",        "%15s",          "fs ref",      "==========", ],
        [   "%10s",      "%10.4f",              "RA",              "==", ],
        [   "%10s",      "%10.4f",             "Dec",             "===", ],
        [    "%6s",         "%6d",            "tmax",            "====", ],
        [    "%6s",         "%6d",           "seeds",           "=====", ],
        [   "%18s", "%9.2e±%8.2e", "flux [ph/cm²/s]", "===============", ],
        [   "%16s", "%9.2f±%6.2f",              "TS",              "==", ],
      # [    "%6s",       "%6.2f",             "√TS",             "===", ],
        [   "%15s", "%8.2f±%6.2f",            "N_on",            "====", ],
        [   "%15s", "%8.2f±%6.2f",           "N_off",           "=====", ],
        [   "%15s", "%8.2f±%6.2f",             "N_s",             "===", ],
        [   "%11s", "%6.2f±%4.2f",           "Li&Ma",           "=====", ],
        [    "%7s",       "%7.4f",           "alpha",           "=====", ],
    ]

    header_fmt = " ".join([r[0] for r in fields]) # headers format
    values_fmt = " ".join([r[1] for r in fields]) # values format
    print(header_fmt % tuple([r[2] for r in fields])) # titles
    print(header_fmt % tuple([r[3] for r in fields])) # sub_titles separator
    for d in sorted(data, key=lambda i: (-1*i["tmax"], i["ra"], i["dec"])):
        n_seeds = len(d["seed_array"])
        flux_m = array_stats(d["flux_array"])
        ts_m   = array_stats(d["ts_array"])
        N_on_m  = array_stats(d["N_on_array"])
        N_off_m = array_stats(d["N_off_array"])
        N_s_m   = array_stats(d["N_s_array"])
        li_ma_m = array_stats(d["li_ma_array"])
        alpha_m = array_stats(d["alpha_array"]) # useless
        print(values_fmt % (d["name"], d["ra"], d["dec"], d["tmax"], n_seeds,
            flux_m["mean"], flux_m["stdev"],
            ts_m["mean"], ts_m["stdev"],
            # sqrt_ts_m["mean"],
            N_on_m["mean"],  N_on_m["stdev"],
            N_off_m["mean"], N_off_m["stdev"],
            N_s_m["mean"],   N_s_m["stdev"],
            li_ma_m["mean"], li_ma_m["stdev"],
            alpha_m["mean"]))

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
    A, mu, sigma, Q = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * 1. / (2. * np.pi * sigma**2)* np.exp(exp_num / exp_den) + Q

def dynamic_bin_number(arr, max_val=None):
    n = max(arr)-min(arr)
    if max_val is not None and n > max_val:
        n = max_val
    return int(n)

# seaborn graph with distplot. Same data, same gaussian loc/scale
#   import seaborn as sns, numpy as np
#   print(array_stats(d["ts_array"]))
#   print(norm.fit(d["ts_array"]))
#   sns.distplot(d["ts_array"], bins=50, kde=True, fit=norm, norm_hist=False)# , norm_hist=True) #array, bins=n_bins, fit=norm, norm_hist=True
def create_hist(ax, data_arr, xlabel=None, color="blue", dyn_bins=False, density_flag=True):
    n_bins = 50
    if dyn_bins:
        n_bins = dynamic_bin_number(data_arr)
    sts = array_stats(data_arr)

    counts_hist, bins_edges, bin_index_not_used = stats.binned_statistic(data_arr, data_arr, statistic='count', bins=n_bins)
    bins_width = np.array(np.diff(bins_edges), float)
    bins_centres = (bins_edges[:-1] + bins_edges[1:])/2

    counts_hist_normalized = counts_hist / bins_width / np.sum(counts_hist)

    fit_params_norm, pvalue_err_norm = fitting_data(gauss, initial_params=[1., sts['mean'], sts['stdev'], 0.], x=bins_centres, y=counts_hist_normalized, verbosity=False, name=xlabel)
    fitted_hist_normalized = gauss(bins_centres, *fit_params_norm)

    if density_flag:
        ax.bar(bins_centres, height=counts_hist_normalized, width=bins_width, alpha=0.5, edgecolor=color, color=color, label='data')
        ax.plot(bins_centres, stats.norm.pdf(bins_centres, sts["mean"], sts["stdev"]), color="orange", linestyle="--", alpha=0.9, label='stats.norm\nμ:{0:.2e}\nσ:{1:.2e}'.format(sts['mean'], sts['stdev']))
        ax.plot(bins_centres, fitted_hist_normalized, label='curve fit\nμ:{0:.2e}\nσ:{1:.2e}'.format(fit_params_norm[1], fit_params_norm[2]), linestyle="-.", color="green")
        ax.legend(loc='best')
    else:
        # ax.bar(bins_centres, height=counts_hist, width=bins_width, alpha=0.5, edgecolor=color, color=color, label='data')
        counts, bins, patches = ax.hist(data_arr, bins=n_bins, alpha=0.5, edgecolor=color, color=color)

    ax.axvline(sts["mean"], color="orange", linestyle="--", alpha=0.9)
    # ax.axvline(sts["median"], color="black",  linestyle="-.", alpha=0.5)

    if xlabel is not None:
        ax.set_xlabel('{} (bins: {})'.format(xlabel, n_bins))

def plot_data_summary(data):
    rows_num=3
    for d in data:
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=2, figsize=(9, rows_num*3))
        fig.suptitle(d["name"]+" "+str(d["tmax"]), va="top", ha="center")

        create_hist(ax[0][0], d["ts_array"],    color="blue",   xlabel="TS")
        create_hist(ax[1][0], d["li_ma_array"], color="cyan",   xlabel="Li&Ma")
        create_hist(ax[2][0], d["flux_array"],  color="green",  xlabel="Flux [ph/cm²/s]" )
        create_hist(ax[0][1], d["N_on_array"],  color="magenta", xlabel="N on",     dyn_bins=True)
        create_hist(ax[1][1], d["N_off_array"], color="red",    xlabel="N off",    dyn_bins=True)
        create_hist(ax[2][1], d["N_s_array"],   color="yellow", xlabel="N signal", dyn_bins=True)

        # Important: first tight_layout(), after adjust for the title
        fig.tight_layout()
        fig.subplots_adjust(top=0.95)
        # plt.show()
        plt.savefig("%s_%04d.png" % (d["name"], d["tmax"]), format="png")
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

    # { 'name': 'dec_0.5', 'seed': '0', 'tmax': '1800',
    #   'ts': '8114.572', 'on_count': '6997.0', 'off_count': '5268.0',
    #   'alpha': '0.25001633167266846', 'li_ma': '90.08091872622624',
    #   'flux': '3.0695648428928983e-09', 'eflux': '4.839076212803621e-10' }
    ds = data_summary(original_data)
    if not data_summary_is_ok(ds):
        exit(1)

    print_data_summary(ds.values())
    plot_data_summary(list(ds.values()))

