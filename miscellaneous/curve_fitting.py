import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit

# Example from:
#    https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.random.normal.html#numpy.random.normal

# chisquare test with a bin with 0 value is no good (variance is zero) so maybe we can use Fischer test
def binning_data(x, my_bins=10):
    counts_hist, bins_edges, not_used = stats.binned_statistic(x, x, statistic='count', bins=my_bins)
    means_hist, not_used, not_used = stats.binned_statistic(x, x, statistic='mean', bins=my_bins)
    std_hist, not_used, not_used = stats.binned_statistic(x, x, statistic='std', bins=my_bins)

    bins_centres = (bins_edges[:-1] + bins_edges[1:])/2
    bins_width = np.array(np.diff(bins_edges), float)
    
    return {
        'bins_edges': bins_edges,
        'bins_centres': bins_centres,
        'bins_width': bins_width,
        'hist_counts': counts_hist,
        'hist_means': means_hist,
        'hist_std': std_hist,
    }

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html#scipy.optimize.curve_fit
def fitting_data(curve_fn, initial_params=[], bins_centres=[], hist=[], more_output=False):
    res = curve_fit(curve_fn, bins_centres, hist, p0=initial_params, full_output=more_output)
    coeff, var_matrix = res[:2]
    if (len(res) > 2):
        infodict, errmsg, ier = res[2:]
        print("infodict: {}\nerrmsg: {}\nier: {}".format(infodict, errmsg, ier))
    perr = np.sqrt(np.diag(var_matrix)) # Standard deviation errors
    print("Curve fit params:")
    print("{0:>10s}  {1:9s}  {2:9s}".format("param no.", "value", "error"))
    for i, c in enumerate(coeff):
        print("{0:10d}  {1:+8.6f}  {2:+8.6f}".format(i, c, perr[i]))
    return coeff, perr

# vignali: (O - E)² / σ²
def chisquare_test(observed, extimated, variance, coeff_num=0):
    bins_num = len(observed)
    if bins_num != len(extimated) or bins_num != len(variance):
        raise Exception("Need same length for observed, extimated and variance array")
    chi2_ts = np.sum((observed - extimated)**2 / variance)
    chi2_pvalue = None
    my_ddof = bins_num - coeff_num
    print("""
        Χ² test: {}
         pvalue: {}
           ddof: {}""".format(chi2_ts, chi2_pvalue, my_ddof))

    return
    # # 
    # ss_res = np.sum((hist - hist_fit)**2)
    # ss_tot = np.sum((hist - np.mean(hist))**2)
    # r2 = 1 - (ss_res / ss_tot)
    # print("R2:", r2)

# removing zero variances
def chisquare_test_fixed(observed, extimated, variance, coeff_num=0):
    bins_num = len(observed)
    if bins_num != len(extimated) or bins_num != len(variance):
        raise Exception("Need same length for observed, extimated and variance array")

    mask = (np.array(variance) > 0.)
    tmp_observed  = np.array(observed)[mask]
    tmp_extimated = np.array(extimated)[mask]
    tmp_variance  = np.array(variance)[mask]
    # print(tmp_observed)
    # print(tmp_extimated)
    # print(np.mean(tmp_observed-tmp_extimated))
    # print(np.mean((tmp_observed-tmp_extimated)**2))
    # print(np.mean(tmp_variance))

    bins_num = len(tmp_observed)
    my_ddof = bins_num - coeff_num
    chi2_ts = np.sum((tmp_observed - tmp_extimated)**2 / tmp_variance)
    chi2_pvalue = None
    print("""fixed
        Χ² test: {}
         pvalue: {}
           ddof: {}""".format(chi2_ts, chi2_pvalue, my_ddof))

# classic: stats.chisquare (O - E)² / E
# def classic_chisquare_test(observed, extimated, coeff_num=0):
#     bins_num = len(observed)
#     if bins_num != len(extimated):
#         raise Exception("Need same length for observed and, extimated array")
# 
#     my_ddof = bins_num - coeff_num
#     chi2_ts = np.sum( (observed-extimated)**2 / extimated )
#     chi2_pvalue = None
#     print("""classic
#         Χ² test: {}
#          pvalue: {}
#            ddof: {}""".format(chi2_ts, chi2_pvalue, my_ddof))
    
# density plot is a smoothed continuous version of a histogram estimated from data.
# common form of estimation is kernel density estimation
# Plot the distribution with a histogram and maximum likelihood gaussian distribution fit:
#   ax = sns.distplot(data, fit=stats.norm, kde=False)
def seaborn_distplot(data, density_fit_values, bins_centres, my_bins=10, name='Unknown'):
    ax = sns.distplot(data, bins=my_bins, fit=stats.norm, fit_kws={'label':'stats.norm', 'color': 'orange', 'linestyle':'--'}, kde=False, hist=True, norm_hist=False, label='data')
    ax.plot(bins_centres, density_fit_values, label='curve fit', linestyle="-.", color="green")

    plt.gca().set(title='{} data dist (Seaborn)'.format(name), ylabel='Normalized freq')
    plt.legend(loc='best')
    plt.show()

def plot_all(x, barplot, fit_values=None, norm_values=None, name=None):
    plt.bar(x, height=barplot['height'], width=[w*0.99 for w in barplot['width']], alpha=0.4, label="data")
    if fit_values is not None:
        plt.plot(x, fit_values, label='Fit w/ curve_fit', linestyle="-.", color="green")
    if norm_values is not None:
        plt.plot(x, norm_values, label='stats.norm.pdf', linestyle='--', color='orange')
    if name:
        plt.gca().set(title='{} data dist'.format(name))
    plt.gca().set(ylabel='Normalized freq')
    plt.legend(loc='best')
    plt.show()

def raw_data_stats(x):
    nobs, (smin, smax), smean, svar, sskewness, skurtosis = stats.describe(x)
    # Delta Degrees of Freedom (divisor is N - ddof)
    # https://stats.stackexchange.com/questions/3931/intuitive-explanation-for-dividing-by-n-1-when-calculating-standard-deviation
    ddof_1_sigma = x.std(ddof=1)
    if np.sqrt(svar) - ddof_1_sigma > 0.00001:
        raise Exception("sigma haven't ddof=1 ({} <> {})".format(np.sqrt(svar), ddof_1_sigma))
    print("""Data stats (with scipy.stats):
          length: {}
       min / max: {} / {}
            mean: {}
             var: {}
           stdev: {} (ddof=1)
       skewnewss: {}
        kurtosis: {}""".format(nobs, smin, smax, smean, svar, np.sqrt(svar), sskewness, skurtosis))
    print("Normal fit from data (mu, sigma w/ ddof=0):", stats.norm.fit(x))

def gauss(x, *p):
    A, mu, sigma = p
    exp_num = -1*(x-mu)**2
    exp_den = 2.*sigma**2
    return A * np.exp(exp_num / exp_den)

def gauss2(x, *p):
    mu, sigma = p
    norm_factor = 1 / np.sqrt(2.* np.pi * sigma**2) 
    exp_num = -1*(x-mu)**2
    exp_den = 2.*sigma**2
    return norm_factor * np.exp(exp_num / exp_den )

def density_normalization(data, bins_width):
    return data / bins_width / np.sum(data)

# this method is useful just ot keep things together
def chisquare_all_analysis(real_data, expected_data, bins_num=0, data_variance=[], params_num=0):
    my_ddof=bins_num-params_num
    scipy_chi2 = stats.chisquare(real_data, f_exp=expected_data, ddof=my_ddof) # Scipy X²
    print("""Scipy
        Χ² test: {}
         pvalue: {}
           ddof: {}""".format(scipy_chi2[0], scipy_chi2[1], my_ddof))

    chisquare_test_fixed(observed=real_data, extimated=expected_data, variance=data_variance, coeff_num=params_num) # fixed X²
    # classic X² =~ scipy X²
    # classic_chisquare_test(observed=real_data, extimated=expected_data, coeff_num=params_num)

# 'bins_edges': bins_edges,
# 'bins_centres': bins_centres,
# 'bins_width': bins_width,
# 'hist_counts': counts_hist,
# 'hist_means': means_hist,
# 'hist_std': std_hist,
if __name__ == '__main__':
    bins=100
    n_samples = 15000
    arr = [
        { 'name': 'Normal',   'data': np.random.normal(size = n_samples)   },
        { 'name': 'Rayleigh', 'data': np.random.rayleigh(size = n_samples) },
        { 'name': 'Random',   'data': np.random.random(size = n_samples)*8-4 },
        # { 'name': 'Exponential', 'data': np.random.exponential(size = n_samples) },
        # { 'name': 'Laplace',  'data': np.random.laplace(size = n_samples) },
    ]
    for d in arr:
        print('===== {} ====='.format(d['name']))
        raw_data_stats(d['data'])
        binned_data = binning_data(d['data'], bins)
        hist_counts = binned_data['hist_counts']

        # searching gaussian curve fit against data count histogram
        print("Fitting curve vs histogram counts")
        fit_params, pvalue_err = fitting_data(gauss, initial_params=[1., 0., 1.], bins_centres=binned_data['bins_centres'], hist=hist_counts, more_output=False)
        hist_fit = gauss(binned_data['bins_centres'], *fit_params)

        chisquare_all_analysis(hist_counts, hist_fit, bins_num=bins, params_num=len(fit_params), data_variance=[ s**2 for s in binned_data['hist_std'] ])

        hist_counts_n = density_normalization(hist_counts, binned_data['bins_width'])
        print("Fitting curve vs normalized histogram counts")
        fit_params_n, pvalue_err_n = fitting_data(gauss, initial_params=[1., 0., 1.], bins_centres=binned_data['bins_centres'], hist=hist_counts_n, more_output=False)
        hist_fit_n = gauss(binned_data['bins_centres'], *fit_params_n)

        chisquare_all_analysis(hist_counts_n, hist_fit_n, bins_num=bins, params_num=len(fit_params), data_variance=[ s**2 for s in binned_data['hist_std'] ])

        norm_mean, norm_std = stats.norm.fit(d['data'])
        norm_y = stats.norm.pdf(binned_data['bins_centres'], norm_mean, norm_std)
        plot_all(x=binned_data['bins_centres'], barplot={'height': hist_counts_n, 'width': binned_data['bins_width'] }, fit_values=hist_fit_n, norm_values=norm_y, name=d['name'])
        # seaborn_distplot(d['data'], density_fit_values=density_hist_fit, bins_centres=binned_data['bins_centres'], my_bins=bins, name=d['name'])

