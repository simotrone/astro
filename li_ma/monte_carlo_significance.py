from scipy import average as avg
from scipy.stats import poisson, norm
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pprint
import lib.significance as li_ma

def describe_data(data, title='No title', **opts):
    fmt  = '      title: "{title}"\n'
    fmt += '     length: {length}\n'
    fmt += '   min, max: {min} - {max}\n'
    fmt += '    average: {avg:.2f}\n'
    if 'mu' in opts:
        fmt += '   expected: {mu} (mu)\n'
    if 'alpha' in opts:
        fmt += '      alpha: {alpha:.2f}\n'
    print(fmt.format(title=title, length=len(data), min=min(data), max=max(data), avg=avg(data), **opts))

def show_random_simulation(data, title=None, mu=None, alpha=None, plot=True):
    describe_data(data, title=title, alpha=alpha, mu=mu)

    if plot:
        fig, ax = plt.subplots(1, 1, tight_layout=False)
        N, bins, patches = ax.hist(data, bins=20, alpha=0.6, rwidth=0.95)
        ax.set_ylabel('counts')
        ax.set_xlabel('random values')
        if title:
            ax.set_title(title)

        text = f'N: {len(data)}'
        if mu:
            text += '\n'+r'$<\mu>$: '+str(mu)
            ax.axvline(mu, ymin=0, ymax=0.95, color='r', linestyle='--')
        if alpha:
            text += '\n'+r'$\alpha$: '+f'{alpha:.3f}'

        tpos = [ ax.get_xlim()[1]*0.8, ax.get_ylim()[1]*0.85 ]
        ax.text(tpos[0], tpos[1], text)
        plt.show()

def describe_signal_significances(signals):
        sign_5  = [i['significance_5']  for i in signals]
        sign_9  = [i['significance_9']  for i in signals]
        sign_17 = [i['significance_17'] for i in signals]
        describe_data(sign_5, title='Significance eq.5')
        describe_data(sign_9, title='Significance eq.9')
        describe_data(sign_17, title='Significance eq.17')


def poisson_rng(mu, attempts):
    r = poisson.rvs(mu, size=attempts)
    return r

def compute_signal_and_significance(observations, alpha):
    """
    observation: random data observation as tuple(n_on, n_off) in iterable
    alpha: the alpha rate t_on/t_off
    """
    results = []
    header=True
    for obs in observations:
        n_on, n_off = obs
        signal = n_on - alpha * n_off

        try:
            significance_5 = li_ma.eq_5(n_on, n_off, alpha)
            significance_9 = li_ma.eq_9(n_on, n_off, alpha)
            significance_17 = li_ma.eq_17(n_on, n_off, alpha)
        except ValueError as e:
            print('Error with these data:', n_on, n_off, alpha)
            raise Exception(e)

        results.append({
            'n_on': n_on,
            'n_off': n_off,
            'signal': signal,
            'significance_5': significance_5,
            'significance_9': significance_9,
            'significance_17': significance_17,
        })

        if OPTS.verbose > 1 and len(results) < 5:
            if header:
                fmt = '{:>4s} {:>4s} | {:>7s} | {:>7s} {:>7s} {:>7s}'
                print(fmt.format('on', 'off', 'signal', 'eq.5', 'eq.9', 'eq.17'))
                print(fmt.format(4*'-', 4*'-', 7*'-', 7*'-', 7*'-', 7*'-'))
                header=False
            fmt = '{:4d} {:4d} | {:+7.2f} | {:+7.2f} {:+7.2f} {:+7.2f}'
            print(fmt.format(n_on, n_off, signal,
                significance_5, significance_9, significance_17))

    return results


def compute_integral_frequency_distribution(data, total_observations, keys, thresholds=[0,1,2,3,4,5]):
    result = {}
    for k in keys:
        result[k] = {                           # for each significance equation:
            'count': np.zeros(len(thresholds)), # counts per significance threshold
            'freq': np.zeros(len(thresholds)),  # frequency of threshold significances
            'negative_significance': 0,         # EXTRA counts for negative significances
        }

    for d in data:
        for k in keys:
            significance_val = d[k]
            # count significance values comparing thresholds
            for i, v in enumerate(thresholds):
                if significance_val >= v:
                    result[k]['count'][i] += 1

            # EXTRA negative significances count
            if significance_val < 0:
                result[k]['negative_significance'] += 1

    # compute frequency
    for k in keys:
        for i, c in enumerate(result[k]['count']):
            result[k]['freq'][i] = c / total_observations

    pp = pprint.PrettyPrinter(indent=2)
    if OPTS.verbose > 1:
        pp.pprint(result)

    # EXTRA check for negative significances
    for k in keys:
        if result[k]['negative_significance'] < 1:
            continue
        pp.pprint(result[k])
        raise Exception(f'WARNING: Unexpected negative significances for {k}')

    return result


def plot_frequency(data, alpha=None, n=None, expected_off=None):
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    for p in data:
        ax.plot(p['x'], p['y'],
            linestyle='', color=p['color'], marker=p['marker'],
            label=p['label'], alpha=0.8)

    ax.plot(p['x'], norm.sf(p['x']), color='black', linestyle='-.', alpha=0.6, label='Gaussian probability')
    ax.set_title(r'$\alpha$'+f'={alpha:.2f}, <N off>={expected_off}, N={n}')
    ax.set_yscale('log')
    ax.set_xlabel('Significance')
    ax.set_ylabel('p (≥ S)')
    ax.legend(loc='best')
    plt.show()


def starts_simulation(expected_N_off, alpha, attempts):
    # generate off counts
    random_N_off = poisson_rng(expected_N_off, attempts)
    # generate on counts (only background as on counts)
    expected_N_on = expected_N_off * alpha
    random_N_on = poisson_rng(expected_N_on, attempts)

    if OPTS.verbose > 0:
        print(r'Simulations data with:')
        print(f'       expected N off: {expected_N_off:4d}')
        print(f'       expected N on : {expected_N_on:6.1f}')
        print(f'                alpha: {alpha:6.2f}')
    if OPTS.verbose > 1:
        plot_flag=True if OPTS.verbose > 2 else False
        print('Details:')
        show_random_simulation(random_N_off, title='Random N off', mu=expected_N_off, alpha=alpha, plot=plot_flag)
        show_random_simulation(random_N_on, title='Random N on', mu=expected_N_on, alpha=alpha, plot=plot_flag)

    #####################
    # create observations
    # no observation with on or off counts <= 0
    simulated_observations = [t for t in zip(random_N_on, random_N_off) if t[0] > 0 and t[1] > 0]
    all_computed_signals = compute_signal_and_significance(simulated_observations, alpha)

    ################
    # filter signals
    # only positive signals contribute to frequency
    positive_signals = [s for s in all_computed_signals if s['signal'] > 0]

    if OPTS.verbose > 0:
        print(f'Simulated observation: {len(simulated_observations)}')
        print(f' all computed signals: {len(all_computed_signals)}')
        print(f'with positive signals: {len(positive_signals)}')

    if OPTS.verbose > 1:
        describe_signal_significances(positive_signals)

    # integral frequency distribution
    # the frequency is computed on simulated observations,
    # but only positive signals are elaborated.
    # significance thresholds: [0, 0.5, 1, ..., 4.5, 5]
    thresholds = np.linspace(0, 5, 11)
    resulting_freq = compute_integral_frequency_distribution(
        positive_signals,
        len(simulated_observations),
        keys=['significance_5', 'significance_9', 'significance_17'],
        thresholds=thresholds,
    )

    ###############
    data_to_plot = [
        { 'x': thresholds,
          'y': resulting_freq['significance_5']['freq'],
          'label': 'eq.5',
          'marker': '+',
          'color': 'blue', },
        { 'x': thresholds,
          'y': resulting_freq['significance_9']['freq'],
          'label': 'eq.9',
          'marker': 'x',
          'color': 'green', },
        { 'x': thresholds,
          'y': resulting_freq['significance_17']['freq'],
          'label': 'eq.17',
          'marker': 'o',
          'color': 'orange', }
    ]
    plot_frequency(data_to_plot, alpha=alpha, expected_off=expected_N_off, n=len(simulated_observations))


def main():
    cases = [
        { 'alpha':  0.1, 'expected_N_off': 140, 'attempts': 1_000_000 },
        { 'alpha':  0.5, 'expected_N_off':  28, 'attempts': 1_000_000 },
        { 'alpha':  1.0, 'expected_N_off':  10, 'attempts': 1_000_000 },
        { 'alpha':  1.0, 'expected_N_off':  14, 'attempts': 1_000_000 },
        { 'alpha':  2.0, 'expected_N_off':  14, 'attempts': 1_000_000 },
        { 'alpha': 10.0, 'expected_N_off':  14, 'attempts': 1_000_000 },
    ]

    for c in cases:
        starts_simulation(c['expected_N_off'], c['alpha'], c['attempts'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Monte Carlo simulations for Li&Ma paper")
    parser.add_argument('-v', '--verbose', action='count', default=0)
    OPTS = parser.parse_args()
    main()
