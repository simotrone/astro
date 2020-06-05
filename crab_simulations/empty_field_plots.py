import argparse
import csv
import logging
import os
import numpy as np
import pprint

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

    return result

def plot_frequency(data, n=None, title=None, save=False):
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    for p in data:
        ax.plot(p['x'], p['y'],
            linestyle='', color=p['color'], marker=p['marker'],
            label=p['label'], alpha=0.8)

    ax.plot(p['x'], norm.sf(p['x']), color='black', linestyle='-.', alpha=0.6, label='Gaussian probability')
    if title:
        ax.set_title(title)
    ax.set_yscale('log')
    ax.set_xlabel('Significance threshold')
    ax.set_ylabel('p (≥ S)')
    ax.legend(loc='best')
    if save:
        plt.savefig(save)
    else:
        plt.show()


def significance_analysis(data, opts):
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
          'label': 'photometric Li&Ma',
          'marker': 'o',
          'color': 'orange', },
        { 'x': thresholds,
          'y': resulting_freq['sqrt_ts']['freq'],
          'label': 'ctools sqrt(TS)',
          'marker': 'v',
          'color': 'blue', },
        # { 'x': thresholds,
        #   'y': resulting_freq['li_ma']['freq'],
        #   'label': 'Li&Ma over ctools data',
        #   'marker': 'x',
        #   'color': 'green', },
    ]
    save_filename = None
    if opts.save:
        save_filename = f'empty_field_{opts.tmax:04d}.png'
    plot_frequency(data_to_plot, n=total, save=save_filename,
        title=f't={opts.tmax} sec Np={len(positive_signals)} N={total}')

def main(opts):
    data = get_raw_data(opts.dir, limit=opts.limit)

    # filter data by tmax (required parameter)
    filtered_data = [d for d in data if d['tmax'] == opts.tmax]
    if opts.check_data:
        check_data(filtered_data)

    significance_analysis(filtered_data, opts)



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
    args = parser.parse_args()
    main(args)

