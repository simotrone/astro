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

if __name__ == "__main__":
    simple()
    triple()
    named_data()
    subplot()
    with_text()
    
    exit(0)
