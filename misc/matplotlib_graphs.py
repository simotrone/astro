import matplotlib.pyplot as plt
import numpy as np

# Esempi di grafici
# https://matplotlib.org/3.1.0/tutorials/introductory/pyplot.html

# simple graph with array data
def simple():
    plt.plot([1, 2, 3, 4], [1, 4, 9, 16])
    plt.axis([0, 6, 0, 20]) # axis range
    plt.ylabel('some numbers')
    plt.show()

# triple graph with numeric series
def triple():
    t = np.arange(0., 5., 0.2) # (min, max, step)
    plt.plot(t, t, 'r--')   # tratteggio
    plt.plot(t, t**2, 'bs') # quadretti
    plt.plot(t, t**3, 'g^') # triangoli
    plt.show()

# graph with named data series
def named_data():
    data = {'a': np.arange(50),
            'color': np.random.randint(0, 50, 50),
            'size': np.random.randn(50)}
    data['b'] = data['a'] + 10 * np.random.randn(50)
    data['size'] = np.abs(data['size']) * 100
    print(data)
    # plt.scatter('a', 'b', data=data) # scatterplot (a,b)
    plt.scatter('a', 'b', data=data, c='color', s='size') # opzioni c: color, s: size
    plt.xlabel('entry a')
    plt.ylabel('entry b')
    plt.show()

# subplots!
def subplot():
    names = ['group_a', 'group_b', 'group_c']
    values = [1, 10, 100]
    plt.figure(figsize=(9, 3)) # width, height in inches.
    
    plt.subplot(131)
    plt.bar(names, values) # histo
    
    plt.subplot(132)
    plt.scatter(names, values) # scatter
    
    plt.subplot(133)
    plt.plot(names, values) # line
    
    plt.suptitle('Categorical Plotting')
    plt.show()

def with_text():
    mu, sigma = 100, 15
    x =  mu + sigma * np.random.randn(10_000)
    # print(x, "\nlen: ", len(x)) # 10_000 records
    n, bins, patches = plt.hist(x, 50, density=1, facecolor='g', alpha=0.75)

    plt.xlabel('Smarts')
    plt.ylabel('Probability')
    plt.title('Histogram of IQ')

    plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.axis([40, 160, 0, 0.03])

    plt.grid(True)
    plt.show()

def power_law(x, k_0, g, E_0):
    return k_0 * ( x / E_0 )**g

def plot_cosmic_rays():

    pivot_energy = 1e9 # GeV
    en = np.linspace(1e8,1e20) # 0.1 GeV (9) -> TeV (12) -> PeV (15) -> 100 EeV (18) 
    if False: # wiebel-sooth (1998)
        prefactor = 3.01
        gamma = -2.68
        plt.plot(en, power_law(en, prefactor, gamma, pivot_energy), label='wiebel-sooth')

    if False: # horandel (2003)
        prefactor = 2.16
        gamma = -2.66
        plt.plot(en, power_law(en, prefactor, gamma, pivot_energy), label='horandel')

    if False: # beringer (2012)
        prefactor = 1.8
        gamma = -2.7
        plt.plot(en, power_law(en, prefactor, gamma, pivot_energy), label='beringer')

    # knee, ankle
    x = [
        np.linspace(1e8, 0.5e16),
        np.linspace(0.5e16, 1e17),
        np.linspace(1e17, 0.5e19),
        np.linspace(0.5e19, 1e20),
    ]
    pivots = [ 1e9, 0.73e10, 1.2e10, 0.46e8 ]
    gammas = [ -2.7, -3.1, -3.2, -2.5]
    if False: # clean plot
        for i, g in enumerate(gammas):
            plt.plot(x[i], power_law(x[i], 3.01, g, pivots[i]))

    if True:
        for i, g in enumerate(gammas):
            plt.plot(x[i], power_law(x[i], 3.01, g, pivots[i])*x[i]**2.6)

    plt.xlabel('Energy (eV)')
    plt.ylabel('')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.show()

    

if __name__ == "__main__":
    # simple()
    # triple()
    # named_data()
    # subplot()
    # with_text()

    # prefactor = 5.7e-16
    # gamma = -2.48
    # pivot_energy = 0.3e6

    # cosmic rays:
    plot_cosmic_rays()
    
    exit(0)
