import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from scipy.stats import chi2
from scipy.optimize import curve_fit

# Example from:
#    https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.random.normal.html#numpy.random.normal

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2.html#scipy.stats.chi2
def show_chi2_things():
    df = 55
    mean, var, skew, kurt = chi2.stats(df, moments="mvsk")
    print("""==== Chi2 ====
         mean: {}
          var: {}
    skewnewss: {}
     kurtosis: {}""".format(mean, var, skew, kurt))

    fig, ax = plt.subplots(1, 1)
    # display the probabily density function (pdf):
    x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
    ax.plot(x, chi2.pdf(x, df), 'r-', lw=5, alpha=0.6, label='chi2 pdf')

    # freeze the distribution and display the frozen pdf:
    rv = chi2(df)
    ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

    # check accuracy of cdf and ppf
    vals = chi2.ppf([0.001, 0.5, 0.999], df)
    print("Chi2 test:", np.allclose([0.001, 0.5, 0.999], chi2.cdf(vals, df)))

    r = chi2.rvs(df, size=5000) # random variates
    ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
    ax.legend(loc='best', frameon=False)
    plt.show()

# normalized normal distribution
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-1*(x-mu)**2/(2.*sigma**2))

# def gauss(x, *p):
#     useless, mu, sigma = p
#     return 1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-1*(x-mu)**2/(2.*sigma**2))

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

# density plot is a smoothed continuous version of a histogram estimated from data.
# common form of estimation is kernel density estimation
def seaborn_distplot(x, my_bins=10, name=None):
    hist, bin_edges = np.histogram(x, density=True, bins=my_bins)
    # center the bins
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    initial_params = [1., 0., 1.]
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html#scipy.optimize.curve_fit
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=initial_params)
    print('Curve fit params (A, mu, sigma):', coeff)
    print('Standard deviation errors:', np.sqrt(np.diag(var_matrix)))
    hist_fit = gauss(bin_centres, *coeff)

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
           sigma: {} (ddof=1)
       skewnewss: {}
        kurtosis: {}""".format(nobs, smin, smax, smean, svar, np.sqrt(svar), sskewness, skurtosis))
    print("Normal fit from data:", norm.fit(x))

if __name__ == '__main__':
    bins=100
    n_samples = 15000
    arr = [
        { 'name': 'Normal',   'data': np.random.normal(size = n_samples)   },
        { 'name': 'Rayleigh', 'data': np.random.rayleigh(size = n_samples) },
        { 'name': 'Laplace',  'data': np.random.laplace(size = n_samples) },
        # { 'name': 'Chisquare', 'data': np.random.chisquare(57, size = n_samples) },
    ]
    for d in arr:
        print('===== {} ====='.format(d['name']))
        data_stats(d['data'])
        seaborn_distplot(d['data'], bins, name=d['name'])
    # show_chi2_things()
    # matplot_all_together(data, bins)

