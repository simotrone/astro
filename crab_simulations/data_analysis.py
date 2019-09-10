import argparse
import os
import sys
import statistics
from lib.exporter.csv import CSVExporter as csvex
import matplotlib.pyplot as plt
from scipy.stats import norm

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
def create_hist(ax, data_arr, xlabel=None, color="blue", dyn_bins=False, exp_fmt=False, density_flag=True):
    n_bins = 50
    if dyn_bins:
        n_bins = dynamic_bin_number(data_arr)
    stats = array_stats(data_arr)

    if density_flag:
        counts, bins, patches = ax.hist(data_arr, bins=n_bins, alpha=0.5, edgecolor=color, color=color, density=True)
        ax.plot(bins, norm.pdf(bins, stats["mean"], stats["stdev"]), color="grey", linestyle="--", alpha=0.9)
    else:
        counts, bins, patches = ax.hist(data_arr, bins=n_bins, alpha=0.5, edgecolor=color, color=color)

    ax.axvline(stats["mean"],   color="grey", linestyle="--", alpha=0.9)
    ax.axvline(stats["median"], color="black",  linestyle="-.", alpha=0.5)

    legend_fmt = 'bins=%d\nmean=%.2f\nmedian=%.2f\nσ=%.2f'
    if exp_fmt:
        legend_fmt = 'bins=%d\nmean=%.2e\nmedian=%.2e\nσ=%.2e'
    ax.text(x=0.65, y=0.70, s=legend_fmt % (n_bins, stats["mean"], stats["median"], stats["stdev"]), transform=ax.transAxes)

    # if int is True:
    #     from matplotlib.ticker import FormatStrFormatter
    #     bin_centers = [ (bins[i]+bins[i+1])*0.5 for i in range(len(bins)-1) ]
    #     ax.set_xticks(bin_centers)
    #     #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    #     ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    if xlabel is not None:
        ax.set_xlabel(xlabel)

def plot_data_summary(data):
    rows_num=3
    for d in data:
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=2, figsize=(9, rows_num*3))
        fig.suptitle(d["name"]+" "+str(d["tmax"]), va="top", ha="center")

        create_hist(ax[0][0], d["ts_array"],    color="blue",   xlabel="TS")
        create_hist(ax[1][0], d["li_ma_array"], color="cyan",   xlabel="Li&Ma")
        create_hist(ax[2][0], d["flux_array"],  color="green",  xlabel="Flux [ph/cm²/s]", exp_fmt=True)
        create_hist(ax[0][1], d["N_on_array"],  color="orange", xlabel="N on",     dyn_bins=True)
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




