import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from scipy.stats import chisquare
from scipy.optimize import curve_fit

# Example from:
#    https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.random.normal.html#numpy.random.normal

# normalized normal distribution
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

# density plot is a smoothed continuous version of a histogram estimated from data.
# common form of estimation is kernel density estimation
def seaborn_distplot(x, my_bins=10, name=None):
    hist, bin_edges = np.histogram(x, density=True, bins=my_bins)
    # center the bins
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    initial_params = [1., 0., 1.]
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html#scipy.optimize.curve_fit
    coeff, var_matrix, infodict, errmsg, ier = curve_fit(gauss, bin_centres, hist, p0=initial_params, full_output=True)
    # print("infodict:",infodict, "\nerrmsg:", errmsg, "\nier:", ier)
    perr = np.sqrt(np.diag(var_matrix)) # Standard deviation errors

    print("Curve fit params:")
    print("{0:>10s}  {1:9s}  {2:9s}".format("param no.", "value", "error"))
    for i, c in enumerate(coeff):
        print("{0:10d}  {1:+8.6f}  {2:+8.6f}".format(i, c, perr[i]))
    hist_fit = gauss(bin_centres, *coeff)
    my_ddof=my_bins-1-len(coeff)
    chi2_ts, chi2_pvalue = chisquare(hist, f_exp=hist_fit, ddof=my_ddof)
    print("""
        Χ² test: {}
         pvalue: {}
           ddof: {}""".format(chi2_ts, chi2_pvalue, my_ddof))
    # 
    ss_res = np.sum((hist - hist_fit)**2)
    ss_tot = np.sum((hist - np.mean(hist))**2)
    r2 = 1 - (ss_res / ss_tot)
    print("R2:", r2)

    # Plot the distribution with a histogram and maximum likelihood gaussian distribution fit:
    # ax = sns.distplot(x, fit=norm, kde=False)
    ax = sns.distplot(x, bins=my_bins, fit_kws={'color': 'orange', 'label':'scipy.norm.pdf', 'linestyle':'--'}, kde=False, hist=True, fit=norm, norm_hist=False, label='data')
    # ax.plot(bin_centres, hist, label='Test data', color='blue')
    ax.plot(bin_centres, hist_fit, label='curve fit', linestyle="-.", color="green")


    my_title = 'Seaborn'
    if name:
        my_title = '{} data dist ({})'.format(name, my_title)
    plt.gca().set(title=my_title, ylabel='Normalized freq')
    plt.legend(loc='best')
    plt.show()

def matplot_all_together(x, my_bins=10):
    hist, bin_edges = np.histogram(x, density=True, bins=my_bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    params = [1., 0., 1.] # initial guess for fitting params (A, mu, sigma)
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html#scipy.optimize.curve_fit
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=params)
    print('Fitted params (A, mu, sigma):', coeff)
    hist_fit = gauss(bin_centres, *coeff)
    plt.hist(x, bins=my_bins, density=True, alpha=0.3, label="data")
    # plt.plot(bin_centres, hist, label='Test data')
    norm_mean, norm_std = norm.fit(x)
    print('Norm fit (mu, sigma w/ ddof=0): {}, {}'.format(norm_mean, norm_std))
    plt.plot(bin_centres, norm.pdf(bin_centres, norm_mean, norm_std), label='scipy.norm.pdf', linestyle='--', color='orange')
    plt.plot(bin_centres, hist_fit, label='Fit w/ curve_fit', linestyle="-.", color="green")
    plt.gca().set(title='Matplotlib', ylabel='Normalized freq')
    plt.legend(loc='best')
    plt.show()

def data_stats(x):
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
    print("Normal fit from data:", norm.fit(x))

if __name__ == '__main__':
    bins=100
    n_samples = 15000
    arr = [
        { 'name': 'Normal',   'data': np.random.normal(size = n_samples)   },
        { 'name': 'Rayleigh', 'data': np.random.rayleigh(size = n_samples) },
        # { 'name': 'Laplace',  'data': np.random.laplace(size = n_samples) },
        # { 'name': 'Chisquare', 'data': np.random.chisquare(57, size = n_samples) },
    ]
    for d in arr:
        print('===== {} ====='.format(d['name']))
        data_stats(d['data'])
        seaborn_distplot(d['data'], bins, name=d['name'])
    # show_chi2_things()
    # matplot_all_together(data, bins)

