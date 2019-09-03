import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from scipy.stats import chi2
from scipy.optimize import curve_fit

# Example from:
#    https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.random.normal.html#numpy.random.normal

def show_hist_matplot_lib(x, my_bins=10):
    mu = x.mean()
    sigma = x.std(ddof=1)
    # Delta Degrees of Freedom (divisor is N - ddof)
    # https://stats.stackexchange.com/questions/3931/intuitive-explanation-for-dividing-by-n-1-when-calculating-standard-deviation
    print(            "Mu:", mu)
    print("Sigma (ddof=1):", sigma) # ddof=1
    count, bins, patches = plt.hist(x, bins=my_bins, density=True)
    plt.gca().set(title='Freq histo (matplotlib)', ylabel='Normalized freq')
    # draw gaussian upon histo
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( -1*(bins-mu)**2 / (2*sigma**2) ), linewidth=2, color='r')
    plt.show()

def show_simple_distplot(x):
    sns.distplot(x, color="orange", hist_kws={'alpha':.6}, kde_kws={'linewidth':2})
    plt.show()

def show_distplot_seaborn(x, my_bins=10):
    sns.distplot(x, bins=my_bins, kde=True, fit=norm, norm_hist=False)
    plt.gca().set(title='Seaborn histo', ylabel='Normalized freq')
    plt.show()

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

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-1*(x-mu)**2/(2.*sigma**2))

def best_gaussian_fit(x, my_bins=10):
    hist, bin_edges = np.histogram(x, density=True, bins=my_bins)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    # print("Edges:", bin_edges)
    # print("Centres:", bin_centres)
    params = [1., 0., 1.] # initial guess for fitting params (A, mu, sigma)
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html#scipy.optimize.curve_fit
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=params)
    # print("Coeff (A, mu sigma):", coeff)
    # print("var_matrix:", var_matrix)
    hist_fit = gauss(bin_centres, *coeff)
    plt.plot(bin_centres, hist, label='Test data')
    plt.plot(bin_centres, hist_fit, label='Fitted data')
    print('Fitted params (A, mu, sigma):', coeff)
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    # data = np.random.normal(size = 5000)
    data = np.random.rayleigh(size = 5000)
    
    # stats
    nobs, (smin, smax), smean, svar, sskewness, skurtosis = stats.describe(data)
    print("""
       length: {}
    min / max: {} / {}
         mean: {}
          var: {}
    skewnewss: {}
     kurtosis: {}""".format(nobs, smin, smax, smean, svar, sskewness, skurtosis))
    
    bins=64
    # show_hist_matplot_lib(data, bins)
    # show_simple_distplot(data)
    show_distplot_seaborn(data, bins)
    # show_chi2_things()
    best_gaussian_fit(data, bins)
    print("Normal fit from data:", norm.fit(data))

