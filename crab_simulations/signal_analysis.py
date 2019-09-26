import argparse
import numpy as np
import os
import statistics
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import curve_fit
from lib.exporter.csv import CSVExporter as csvex
import logging

# Example:
# PYTHONPATH=../path/to/lib/ python signal_analysis.py *tsv
# PYTHONPATH=../astro_code/ python signal_analysis.py dec_*/*tsv ra_*/*tsv

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

# { 'name': 'dec_0.5', 'seed': '0', 'tmax': '1800', 'ts': '8114.572',
#   index_value: xxx      index_error: xxx
#   prefactor_value: xxx  prefactor_error: xxx
#   pivot_value: xxx      pivot_error: xxx
#   'flux': '3.0695648428928983e-09', 'eflux': '4.839076212803621e-10'
#   'on_count': '6997.0', 'off_count': '5268.0',
#   excess_count: xxx
#   'alpha': '0.25001633167266846', 'li_ma': '90.08091872622624',
# TODO: è corretto fare la media di li & ma? oppure è meglio calcolarlo sulla
#       media di N_on e N_off
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
                    # index_value
                    # index_error
                    # prefactor_value
                    # prefactor_error
                    # pivot_value
                    # pivot_error
                    'flux':  [],
                    'eflux': [],
                    'N_on':  [],
                    'N_off': [],
                    'N_exc': [],
                    'alpha': [],
                    'li_ma': [],
                },
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
            logging.warning("{0:15s} seed:{1:3d} tmax:{2:4d} ({3:.0f} on, {4:2.0f} off): Negative ts {5:.2f}".format(i["name"], int(i["seed"]), int(i["tmax"]), float(i["on_count"]), float(i["off_count"]), float(i["ts"])))
        elif i["li_ma"] is None:
            logging.warning("{0:15s} seed:{1:3d} tmax:{2:4d} ({3:.0f} on, {4:2.0f} off): Cannot calculate Li&Ma".format(i["name"], int(i["seed"]), int(i["tmax"]), float(i["on_count"]), float(i["off_count"])))
    return o

# WARNING: this function augment the input data struct
def data_augmentation(data, bins_number=50):
    fields = [
        { 'name': 'N_on',  'dyn_bins': True },
        { 'name': 'N_off', 'dyn_bins': True },
        { 'name': 'N_exc', 'dyn_bins': True },
        { 'name': 'li_ma', 'dyn_bins': False },
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
                'bins_edges':   bins_edges,
                'bins_centres': bins_centres,
                'bins_width':   bins_width,
                'fit_coeff':  fit_coeff,
                'pvalue_err': pvalue_err,
            }
    return data

def array_stats(arr):
    stat = {
        'n': len(arr),
        'mean':   statistics.mean(arr),
        'stdev':  statistics.pstdev(arr),
        'median': statistics.median(arr),
    }
    return stat

def print_txt_data_summary(data):
    fields = [
      #   h_format,         v_format,     title,         sub_t
        [   '%15s',           '%15s',  'fs ref',  '==========', ],
        [   '%10s',           '%10s',      'RA',          '==', ],
        [   '%10s',           '%10s',     'Dec',         '===', ],
        [    '%6s',            '%6d',    'tmax',        '====', ],
        [    '%6s',            '%6d',   'seeds',       '=====', ],
        [   '%16s',    '%9.2f±%6.2f',      'TS',          '==', ],
        [   '%15s',    '%8.2f±%6.2f',    'N_on',        '====', ],
        [   '%15s',    '%8.2f±%6.2f',   'N_off',       '=====', ],
        [   '%15s',    '%8.2f±%6.2f',     'N_s',         '===', ],
        [   '%11s',    '%6.2f±%4.2f',   'Li&Ma',       '=====', ],
        [    '%7s',          '%7.4f',   'alpha',       '=====', ],
        [   '%26s', '%10.2f %7.2f %6.2f', 'N_on fitting (A, μ, σ)', '=======', ],
        [   '%23s', '%10.2f %5.2f %5.2f',  'N_on pvalue (A, μ, σ)', '=======', ],
    ]

    header_fmt = ' '.join([r[0] for r in fields]) # headers format
    values_fmt = ' '.join([r[1] for r in fields]) # values format
    print(header_fmt % tuple([r[2] for r in fields])) # titles
    print(header_fmt % tuple([r[3] for r in fields])) # sub_titles separator
    for d in sorted(data.values(), key=lambda i: (-1*i['tmax'], i['ra'], i['dec'])):
        n_seeds = len(d['data']['seed'])
        ts_m    = array_stats(d['data']['ts'])
        N_on_m  = d['stats']['N_on']
        N_off_m = d['stats']['N_off']
        N_exc_m = d['stats']['N_exc']
        li_ma_m = d['stats']['li_ma']
        alpha_m = array_stats(d['data']['alpha']) # useless
        if alpha_m['stdev'] > 0.000001:
            logging.error('Just a check. alpha stdev must be 0. alpha={}'.format(alpha_m))
            exit(1)
        print(values_fmt % (d['name'], d['ra'], d['dec'], d['tmax'], n_seeds,
                            ts_m['mean'],    ts_m['stdev'],
                            N_on_m['mean'],  N_on_m['stdev'],
                            N_off_m['mean'], N_off_m['stdev'],
                            N_exc_m['mean'], N_exc_m['stdev'],
                            li_ma_m['mean'], li_ma_m['stdev'],
                            alpha_m['mean'],
                            d['hist']['N_on']['fit_coeff'][0],
                            d['hist']['N_on']['fit_coeff'][1],
                            abs(d['hist']['N_on']['fit_coeff'][2]),
                            d['hist']['N_on']['pvalue_err'][0],
                            d['hist']['N_on']['pvalue_err'][1],
                            d['hist']['N_on']['pvalue_err'][2] ))

def print_html_data_summary(data):
    from  jinja2 import Template
    t = Template("""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8" />
    <style>
    .alnright { text-align: right; }
    .spaced { float: right; width: 3em; }
    table { width: auto; }
    th, td { padding: 10px; }
    th { background-color: #4CAF50; color: white; }
    th { border-bottom: 1px solid #ddd; }
    tr:nth-child(even) { background-color: #f6f6f6; }
    tr:hover {background-color: #f2f2f2;}
    p { margin: 0; padding 3px 5px; }
    ul.gallery li {
        list-style-type: none;
        float: left;
        border: 1px solid #a2a2a2;
        padding: 1em;
        margin: 1em;
    }
    </style>
    <title>Crab signal</title>
</head>
<body>
    <h2>Crab excess counts with on/off analysis</h2>
    <table id="main_table">
        {% set doing = {} %}
        {% for d in rows %}
            {% if d['name'] not in doing %}
                <tr>
                    <th>name</th>
                    <th>tmax <br/> [sec]</th>
                    <th>Total <br/> seeds</th>
                    <th>On source <br/> counts [ph]</th>
                    <th>Off source <br/> counts [ph]</th>
                    <th>Excess <br/> counts [ph]</th>
                    <th>Li &amp; Ma <br/> significance</th>
                </tr>
            {% endif %}
            {% if doing.update({ d['name']: True }) %} {% endif %}
        <tr>
            <td id="{{ d['name'] }}">{{ d['name'] }} (<a href="#{{ d['img'] }}">plot</a>)</td>
            <td class="alnright">{{ d['tmax'] }}</td>
            <td class="alnright">{{ d['data']['seed']|length }}</td>
            <td class="alnright">
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['stats']['N_on']['mean'], d['stats']['N_on']['stdev']) }}                   <span class="spaced">data</span></p>
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['hist']['N_on']['fit_coeff'][1],  d['hist']['N_on']['fit_coeff'][2]|abs) }} <span class="spaced">fit </span></p>
            </td>
            <td class="alnright">
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['stats']['N_off']['mean'], d['stats']['N_off']['stdev']) }}                  <span class="spaced">data</span></p>
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['hist']['N_off']['fit_coeff'][1], d['hist']['N_off']['fit_coeff'][2]|abs) }} <span class="spaced">fit </span></p>
            </td>
            <td class="alnright">
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['stats']['N_exc']['mean'], d['stats']['N_exc']['stdev']) }}                  <span class="spaced">data</span></p>
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['hist']['N_exc']['fit_coeff'][1], d['hist']['N_exc']['fit_coeff'][2]|abs) }} <span class="spaced">fit </span></p>
            </td>
            <td class="alnright">
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['stats']['li_ma']['mean'], d['stats']['li_ma']['stdev']) }}                  <span class="spaced">data</span> </p>
                <p>{{ '{0:.3f} ± {1:.3f}'.format(d['hist']['li_ma']['fit_coeff'][1], d['hist']['li_ma']['fit_coeff'][2]|abs) }} <span class="spaced">fit </span> </p>
            </td>
        </tr>
        {% endfor %}
    </table>

    <h3>Plots</h3>

    <ul class="gallery">
        {% for d in rows %}
        <li id="{{ d['img'] }}"><img src="{{ d['img'] }}" /> <a href="#main_table">back</a></li>
        {% endfor %}
    </ul>
</body>
</html>
""")
    html = t.render(rows=data.values())
    print(html)
    

def fitting_data(curve_fn, initial_params=[], x=[], y=[], verbosity=False, name=None):
    res = curve_fit(curve_fn, x, y, p0=initial_params, full_output=verbosity)
    coeff, var_matrix = res[:2]
    if (len(res) > 2):
        infodict, errmsg, ier = res[2:]
        print('infodict: {}\nerrmsg: {}\nier: {}'.format(infodict, errmsg, ier))
    perr = np.sqrt(np.diag(var_matrix))
    logging.debug('Curve fit params: {}'.format(name))
    logging.debug('{0:>10s}  {1:9s}  {2:9s}'.format('param no.', 'value', 'error'))
    for i, c in enumerate(coeff):
        logging.debug('{0:10d}  {1:+8.6e}  {2:+8.6e}'.format(i, c, perr[i]))
    return coeff, perr

# no good with the quote
def gauss(x, *params):
    A, mu, sigma = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * np.exp(exp_num / exp_den)

def gauss2(x, *params):
    A, mu, sigma = params
    exp_num = -1 * (x-mu)**2
    exp_den = 2. * sigma**2
    return A * 1. / (2. * np.pi * sigma**2)* np.exp(exp_num / exp_den)

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

    # histogram
    ax.bar(bins_centres, height=counts_hist, width=bins_width, alpha=0.5, edgecolor=color, color=color, label='data')

    # normal stats
    # ax.plot(bins_centres, stats.norm.pdf(bins_centres, data_stats["mean"], data_stats["stdev"]), color="orange", linestyle="--", alpha=0.9, label='stats.norm\nμ:{0:.2e}\nσ:{1:.2e}'.format(data_stats['mean'], data_stats['stdev']))
    # ax.axvline(data_stats["mean"], color="blue", linestyle="--", alpha=0.9)

    # gauss fit
    ax.plot(bins_centres, fitted_hist, linestyle="-.", color="green")

    ax.legend( [ Patch(facecolor=color, edgecolor=color),
                 Line2D([0],[0], color="green", linestyle="-.") ],
               [ 'Data. bins:{0}\nμ:{1:.2f}\nσ:{2:.2f}'.format(data['n_bins'], data_stats['mean'], data_stats['stdev']),
                 'Curve fit\nμ:{0:.2f}\nσ:{1:.2f}'.format(fit_params[1], abs(fit_params[2])) ],
               loc='best' )

    ax.set_xlabel('{0}'.format(xlabel))
    ax.set_ylabel('counts')

def plot_data_summary(data, save_img=False):
    rows_num=2
    cols_num=2
    img_format = 'png'
    for d in list(data.values()):
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=cols_num, figsize=(cols_num*4.5, rows_num*4.5))
        fig.suptitle('{} t_{} sec'.format(d["name"], str(d["tmax"])), va="top", ha="center")

        create_hist(ax[0][0], d['hist']['N_on'],  d['stats']['N_on'],  color="magenta", xlabel="N on")
        create_hist(ax[0][1], d['hist']['N_off'], d['stats']['N_off'], color="red",     xlabel="N off")
        create_hist(ax[1][0], d['hist']['N_exc'], d['stats']['N_exc'], color="yellow",  xlabel="N excess")
        create_hist(ax[1][1], d['hist']['li_ma'], d['stats']['li_ma'], color="orange",  xlabel="Li & Ma significance")

        # Important: first tight_layout(), after adjust for the title
        fig.tight_layout()
        fig.subplots_adjust(top=0.90)
        if save_img:
            img_dir = "imgs_signal"
            img_filename = "{0}/signal_{1}_{2:04d}.{3}".format(img_dir, d["name"], d["tmax"], img_format)
            try:
                os.makedirs(img_dir)
            except FileExistsError as e:
                logging.debug("The imgs dir {} already exists".format(img_dir))
            plt.savefig(img_filename, format=img_format)
            plt.close()
            logging.debug("saving {}".format(img_filename))
            d['img'] = img_filename
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

    logging.warning("File read: {}".format(file_count))
    logging.warning("     data: {}".format(len(original_data)))

    ds = data_summary(original_data)
    if not data_summary_is_ok(ds, pointings=5, time_slots=8, different_seeds=5000):
        exit(1)
    
    # check the min/max TS
    if False is True:
        for d in ds.values():
            ns = d['data']['N_exc']
            print('{0:10.5f} {1:10.5f}'.format(np.min(ns), np.max(ns)))
        exit(1)

    # inplace augment
    data_augmentation(ds, bins_number=10)
    plot_data_summary(ds, save_img=True)
    # print_txt_data_summary(ds)
    print_html_data_summary(ds)

