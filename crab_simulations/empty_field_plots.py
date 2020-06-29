import argparse
import csv
import logging
import os
import numpy as np
import pprint
import matplotlib.pyplot as plt
from scipy.stats import norm, chi2, binned_statistic

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('sa') # significances analyzer
pp = pprint.PrettyPrinter(indent=2)

def get_files_from_dir(directory):
    ret = []
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        if not os.path.isfile(filepath):
            log.warn(f'object {filepath} is not a file. skip')
            continue
        if not filepath.endswith('.tsv'):
            log.warn(f'file {filepath} is not a tsv file. skip')
            continue
        ret.append(filepath)
    return ret

def get_files_from_dirs(directories):
    ret = []
    for d in directories:
        try:
            filenames = get_files_from_dir(d)
            ret += filenames
        except NotADirectoryError as e:
            log.warn(e)
            log.warn('skip')
    return ret

def headers_check(headers, fields, filepath):
    for h in headers:
        if h not in fields:
            raise Exception(f"header '{h}' in file '{filepath}' is not valid.")

def get_row_data(row, headers, types, filepath):
    data = {}
    for i, h in enumerate(headers):
        val = row[i]
        if val == '' and types[h] is float:
            val = float('nan')
        try:
            data[h] = types[h](val)
        except ValueError as e:
            log.error(f"file '{filepath}', header '{h}': value '{val}'")
            raise e
    return data

def read_csv(filepath, newline='\n', separator='\t', fields=None):
    data = []
    with open(filepath, newline=newline) as fh:
        reader = csv.reader(fh, delimiter=separator)
        headers_row = None
        if fields:
            headers_row = next(reader, None)
            headers_check(headers_row, fields, filepath)
        for row in reader:
            if headers_row:
                d = get_row_data(row, headers_row, fields, filepath)
            else:
                d = row
            data.append(d)
    return data

def get_raw_data(directories, limit=None):
    data_files = get_files_from_dirs(directories)
    log.info(f'files sourced: {len(data_files)}')
    if limit:
        data_files = data_files[0:limit]

    field_types = { 'name': str,
                    'ra': float, 'dec': float,
                    'seed': int, 'tmax': int,
                    'ts': float,
                    'index_value': float, 'index_error': float,
                    'prefactor_value': float, 'prefactor_error': float,
                    'pivot_value': float, 'pivot_error': float,
                    'flux': float, 'eflux': float,
                    'on_count': float, 'off_count': float, 'alpha': float,
                    'excess_count': float, 'li_ma': float,
                    'phm_on': float, 'phm_off': float,
                    'phm_excess': float, 'phm_li_ma': float, }
    ret = []
    for filepath in data_files:
        data = read_csv(filepath, fields=field_types)
        ret += data
    log.debug(f'data extracted: {len(ret)}')
    return ret

def check_data(data):
    count_fields = ['on_count', 'off_count', 'phm_on', 'phm_off' ]
    min_max_results = {
        'ts':        { 'min': 1e10, 'max': -1e10 },
        'li_ma':     { 'min': 1e10, 'max': -1e10 },
        'phm_li_ma': { 'min': 1e10, 'max': -1e10 },
    }

    for d in data:
        for f in count_fields:
            if d[f] >= 0:
                continue
            log.info(f"field '{f}' value '{d[f]}' [{d['seed']}, {d['tmax']}]")
        for f in min_max_results.keys():
            if d[f] < min_max_results[f]['min']:
                min_max_results[f]['min'] = d[f]
            if d[f] > min_max_results[f]['max']:
                min_max_results[f]['max'] = d[f]
    # recap
    for k, v in min_max_results.items():
        log.info(f"'{k:>10s}': {v['min']:+8.3f}, {v['max']:+8.3f}")


def compute_integral_frequency_distribution(data, keys, thresholds=[0,1,2,3,4,5], total=None):
    result = {}
    for k in keys:
        result[k] = {
            'thresholds': thresholds,
            'count': np.zeros(len(thresholds)),
            'freq': np.zeros(len(thresholds)),
            'freq_err': np.zeros(len(thresholds)),
        }
    for d in data:
        for k in keys:
            for i, t_val in enumerate(thresholds):
                if d[k] >= t_val:
                    result[k]['count'][i] += 1

    # if total => freq
    if total is not None and total > 0:
        for k in keys:
            for i, c in enumerate(result[k]['count']):
                result[k]['freq'][i] = c / total
                result[k]['freq_err'][i] = np.sqrt(c) / total 

    return result

def plot_frequency(data, n=None, title=None, save=False):
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    for p in data:
        ax.errorbar(x=p['x'], y=p['y'],
            yerr=p['yerr'],
            linestyle='', color=p['color'], marker=p['marker'],
            label=p['label'], alpha=0.8)

    ax.plot(p['x'], norm.sf(p['x']), color='black', linestyle='-.', alpha=0.6, label='Gaussian probability')
    # defined p-value: 0.05
    if True:
        norm_isf = norm.isf(0.05)
        ax.plot([0, norm_isf], [0.05, 0.05], color='red', ls=':', alpha=0.9,
            label=f'p-value=0.05, S={norm_isf:.2f}')
        ax.plot([norm_isf, norm_isf], [0, norm.sf(norm_isf)], color='red',
            ls=':', alpha=0.9)
    if title:
        ax.set_title(title)
    ax.set_yscale('log')
    ax.set_xlabel('Significance threshold')
    ax.set_ylabel('p (≥ S)')
    ax.legend(loc='best')
    if save:
        plt.savefig(save)
        log.info(f'significance p-value plot saved in {save}')
    else:
        plt.show()

def p_value_analysis(data, opts):
    """
    p-value analysis with cumulative freq.
    Only positive signals.
    Negative ctools ts => 0
    """
    # photometrics analysis
    total = len(data)
    positive_signals = [d for d in data if d['phm_excess'] > 0]
    log.info(f'      considered data: {len(data)}')
    log.info(f'with positive signals: {len(positive_signals)}')

    # data augmentation. add sqrt(ts) to data
    for p in positive_signals:
        # if False and p['ts'] < -1:
        #     raise Exception(f"negative ts {p['ts']:.2e} is not so tiny. {p}")
        p['sqrt_ts'] = np.sqrt(p['ts']) if p['ts'] > 0 else 0
    # ########

    thresholds = np.linspace(0, 5, 11)
    resulting_freq = compute_integral_frequency_distribution(
        positive_signals,
        keys=['phm_li_ma', 'li_ma', 'sqrt_ts'],
        thresholds=thresholds,
        total=total,
    )

    data_to_plot = [
        { 'x': thresholds,
          'y': resulting_freq['phm_li_ma']['freq'],
          'yerr': resulting_freq['phm_li_ma']['freq_err'],
          'label': 'photometric Li&Ma',
          'marker': '.',
          'color': 'orange', },
        { 'x': thresholds,
          'y': resulting_freq['sqrt_ts']['freq'],
          'yerr': resulting_freq['sqrt_ts']['freq_err'],
          'label': 'ctools sqrt(TS)',
          'marker': '.',
          'color': 'blue', },
        # { 'x': thresholds,
        #   'y': resulting_freq['li_ma']['freq'],
        #   'label': 'Li&Ma over ctools data',
        #   'marker': 'x',
        #   'color': 'green', },
    ]
    save_filename = None
    if opts.save:
        save_filename = f'empty_field_norm_{opts.tmax:04d}.png'
    plot_frequency(data_to_plot, n=total, save=save_filename,
        title=f't={opts.tmax} sec Np={len(positive_signals)} N={total}')

def p_value_analysis_fulldata(data, opts):
    """
    p-value analysis with cumulative freq.
    All data are considered, but we filter with ctools ts >= 0.
    """
    total = len(data)
    accepted_signals = [d for d in data if d['phm_excess'] > 0 or d['ts'] >= 0]
    log.info(f'      considered data: {len(data)}')
    log.info(f'with positive signals: {len(accepted_signals)}')

    # data augmentation. add sqrt(ts) to data
    for p in accepted_signals:
        p['sqrt_ts'] = np.sqrt(p['ts'])

    thresholds = np.linspace(0, 5, 11)
    resulting_freq = compute_integral_frequency_distribution(
        accepted_signals,
        keys=['phm_li_ma', 'sqrt_ts'],
        thresholds=thresholds,
        total=total,
    )

    data_to_plot = [
        { 'x': thresholds,
          'y': resulting_freq['phm_li_ma']['freq'],
          'label': 'photometric Li&Ma',
          'marker': 'o',
          'color': 'orange', },
        { 'x': thresholds,
          'y': resulting_freq['sqrt_ts']['freq'],
          'label': 'ctools sqrt(TS)',
          'marker': 'v',
          'color': 'blue', },
    ]
    save_filename = None
    if opts.save:
        save_filename = f'empty_field_fulldata_norm_{opts.tmax:04d}.png'
    plot_frequency(data_to_plot, n=total, save=save_filename,
        title=f't={opts.tmax} sec Np={len(accepted_signals)} N={total}')

def data_binning(values, bins=10):
    counts, bin_edges, bin_index_not_used = binned_statistic(values, values, statistic='count', bins=bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    bin_widths = np.diff(bin_edges)
    total = len(values)
    freq = counts/total

    ret = {
        'total': total,
        'counts': counts,
        'freq': freq,
        'bin_edges': bin_edges,
        'bin_centres': bin_centres,
        'bin_widths': bin_widths,
        'n_bins': len(bin_edges)-1,
        'max': np.max(values),
        'min': np.min(values),
    }
    # pp.pprint(ret)
    return ret

def plot_counts(data, xlabel=None, title=None, save=False):
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    for p in data:
        ax.errorbar(x=p['x'], y=p['y'],
            xerr=p['xerr'], yerr=p['yerr'],
            linestyle='', color=p['color'], marker=p['marker'],
            label=p['label'], alpha=0.8)
    #ax.errorbar(x=data['bin_centres'], y=data['counts'], fmt='+')
    ax.plot(p['x'], chi2.pdf(p['x'], df=1),
        label='Χ²', linestyle='-.', color='black', alpha=0.6)
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.set_yscale('log')
    ax.set_ylabel('counts (normalized)')
    ax.legend(loc='best')
    if save:
        plt.savefig(save)
        log.info(f'ts distribution plot saved in {save}')
    else:
        plt.show()

def ts_distribution(data, opts):
    all_ts = []
    all_phm_li_ma2 = []
    for d in data:
        # if d['phm_excess'] <= 0: continue
        ts_val = d['ts'] if d['ts'] > 0 else 0
        all_ts.append(ts_val)
        all_phm_li_ma2.append(d['phm_li_ma']**2)

    bins_edges=np.linspace(0, 30, 31)
    ts_binned = data_binning(all_ts, bins=bins_edges)
    phm_binned = data_binning(all_phm_li_ma2, bins=bins_edges)

    data_to_plot = [
        { 'x': phm_binned['bin_centres'],
          'y': phm_binned['freq'],
          'xerr': phm_binned['bin_widths']/2,
          'yerr': np.sqrt(phm_binned['freq']/phm_binned['total']),
          'label': 'Li&Ma²',
          'marker': '.',
          'color': 'orange', },
        { 'x': ts_binned['bin_centres'],
          'y': ts_binned['freq'],
          'xerr': ts_binned['bin_widths']/2,
          'yerr': np.sqrt(ts_binned['freq']/ts_binned['total']),
          'label': 'ctools TS',
          'marker': '.',
          'color': 'blue', },
    ]
    save_filename = None
    if opts.save:
        save_filename = f'empty_field_ts_dist_{opts.tmax:04d}.png'
    plot_counts(data_to_plot, xlabel='Test Statistic', save=save_filename,
        title=f't={opts.tmax} sec N={len(all_ts)}')

def significance_distribution(data, opts):
    all_ts = []
    all_phm_li_ma = []
    for d in data:
        # if d['phm_excess'] <= 0: continue
        all_ts.append(d['ts'])
        all_phm_li_ma.append(d['phm_li_ma'])
    ts_binned = data_binning(all_ts, bins=30)
    phm_binned = data_binning(all_phm_li_ma, bins=30)

    factor=0.7
    nrows=2
    ncols=2
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True,
        figsize=(ncols*6.4*factor, nrows*4.8*factor))

    axes[0][0].errorbar(x=ts_binned['bin_centres'], y=ts_binned['counts'],
        xerr=ts_binned['bin_widths']/2, yerr=np.sqrt(ts_binned['counts']),
        linestyle='', color='blue', marker='.', label='ctools TS')
    axes[0][0].set_title(f't={opts.tmax} sec N={len(all_ts)}')
    axes[0][0].set_xlabel('ctools TS')
    axes[0][0].set_ylabel('counts')

    axes[0][1].errorbar(x=phm_binned['bin_centres'], y=phm_binned['counts'],
        xerr=phm_binned['bin_widths']/2, yerr=np.sqrt(phm_binned['counts']),
        linestyle='', color='orange', marker='.', label='Li&Ma')
    axes[0][1].set_title(f't={opts.tmax} sec N={len(all_phm_li_ma)}')
    axes[0][1].set_xlabel('Significance (phm Li&Ma)')
    axes[0][1].set_ylabel('counts')

    axes[1][0].errorbar(x=ts_binned['bin_centres'], y=ts_binned['freq'],
        xerr=ts_binned['bin_widths']/2,
        yerr=np.sqrt(ts_binned['freq']/ts_binned['total']),
        linestyle='', color='blue', marker='.', label='ctools TS')
    axes[1][0].set_title(f't={opts.tmax} sec N={len(all_ts)}')
    axes[1][0].set_xlabel('ctools TS')
    axes[1][0].set_ylabel('counts (normalized)')
    axes[1][0].set_yscale('log')
    axes[1][0].plot(ts_binned['bin_centres'],
        chi2.pdf(ts_binned['bin_centres'], df=1), label='Χ² pdf',
        linestyle='-.', color='black', alpha=0.6)

    axes[1][1].errorbar(x=phm_binned['bin_centres'], y=phm_binned['freq'],
        xerr=phm_binned['bin_widths']/2,
        yerr=np.sqrt(phm_binned['freq']/phm_binned['total']),
        linestyle='', color='orange', marker='.', label='Li&Ma')
    axes[1][1].set_title(f't={opts.tmax} sec N={len(all_phm_li_ma)}')
    axes[1][1].set_xlabel('Significance (phm Li&Ma)')
    axes[1][1].set_ylabel('counts (normalized)')
    axes[1][1].plot(phm_binned['bin_centres'],
        norm.pdf(phm_binned['bin_centres']), label='norm pdf',
        linestyle='-.', color='black', alpha=0.6)

    for ax in axes.flatten():
        ax.legend(loc='best')

    if opts.save:
        save_filename = f'empty_field_s_dist_{opts.tmax:04d}.png'
        plt.savefig(save_filename)
        log.info(f'significance distribution plot saved in {save_filename}')
    else:
        plt.show()

def export_data(data, fields):
    headers = ['seed'] + fields
    print(",".join(headers))
    for d in data:
        row = [str(d['seed'])]
        for f in fields:
            val = d[f]
            if not isinstance(val, float):
                raise Exception(f"Issue with {val}, it is not float.")
            row.append(str(val))
        print(",".join(row))

def main(opts):
    data = get_raw_data(opts.dir, limit=opts.limit)

    # filter data by tmax (required parameter)
    filtered_data = [d for d in data if d['tmax'] == opts.tmax]
    if opts.check_data:
        check_data(filtered_data)

    if opts.extract:
        fields = ['ts', 'phm_on', 'phm_off', 'phm_excess', 'phm_li_ma']
        export_data(filtered_data, fields=fields)
        exit(0)

    # show distribution of raw significance and TS
    significance_distribution(filtered_data, opts)
    # show distribution of ctools TS and Li&Ma²
    ts_distribution(filtered_data, opts)
    # show cumulative freq vs survival
    p_value_analysis(filtered_data, opts)
    # p_value_analysis_fulldata(filtered_data, opts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze significance from tsv")
    parser.add_argument('dir', help='the dir with tsv files', nargs='+')
    parser.add_argument('-tmax', '--tmax', type=int, required=True,
        help='select data for specific tmax simulation')
    parser.add_argument('-l', '--limit', default=None, type=int,
        help='the number of records to analyze')
    parser.add_argument('--save', action='store_true',
        help='save the plot')
    parser.add_argument('--check-data', action='store_true',
        help='check the raw data')
    parser.add_argument('--extract', action='store_true',
        help='extract data')
    args = parser.parse_args()
    main(args)

