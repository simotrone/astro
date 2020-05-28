from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-5, 5, 100)
zoomed_x = np.linspace(0, 5, 100)

data = [
    { 'x': x, 'fn': norm.pdf, 'c': 'red',  'ls': '-', 'label': 'pdf', 'title': 'Probability density function (pdf)' },
    { 'x': x, 'fn': norm.cdf, 'c': 'blue', 'ls': ':', 'label': 'cdf', 'title': 'Cumulative distribution function (cdf)' },
    { 'x': x, 'fn': norm.sf, 'c': 'black', 'ls': '-.', 'label': 'sf', 'title': 'Survival function (1-cdf)' },
 #   { 'x': x, 'fn': norm.sf, 'c': 'darkgreen', 'ls': '-.', 'label': 'log sf', 'title': 'Survival function (log)', 'yscale': 'log' },
    { 'x': zoomed_x, 'fn': norm.sf, 'c': 'black', 'ls': '-.', 'label': 'survival', 'title': 'Survival zoom [0,5]', 'yscale': 'log' },
]

# plots
factor = 1.7
rows,cols=2,2
fig, axes = plt.subplots(rows, cols,
    figsize=(cols*3*factor, rows*2*factor),
    tight_layout=False)
fig.suptitle('Gaussian things')

ax = axes.flatten()

for i,d in enumerate(data):
    ax[i].plot(d['x'], d['fn'](d['x']), color=d['c'], linestyle=d['ls'], label=d['label'])
    ax[i].set_title(d['title'])
    if 'yscale' in d:
        ax[i].set_yscale(d['yscale'])

# sigma on pdf
for i in range(1,6):
    for j in -i, i:
        ax[0].plot((j,j), (0, norm.pdf(j)), color='gray', linestyle='-', alpha=0.5)

    ax[0].text(x=2, y=norm.pdf(0)*0.6,
        s='1σ: {:.4%}\n'.format(norm.cdf(1)-norm.cdf(-1)) \
         +'2σ: {:.4%}\n'.format(norm.cdf(2)-norm.cdf(-2)) \
         +'3σ: {:.4%}\n'.format(norm.cdf(3)-norm.cdf(-3)) \
         +'4σ: {:.4%}\n'.format(norm.cdf(4)-norm.cdf(-4)) \
         +'5σ: {:.4%}'.format(norm.cdf(5)-norm.cdf(-5)))

# annotation on cumulative
for i in range(0,6):
    if i < 3:
        text_pos = (-5, 0.4+i*0.25)
    else:
        text_pos = ((i-3)*1.4, 0.3+(i-3)*-0.12)
    fmt='{:.'+f'{i}'+'%}'
    ax[1].annotate(fmt.format(norm.cdf(i)),
        xy=(i, norm.cdf(i)),
        xytext=text_pos,
        arrowprops={ 'arrowstyle': '->' },
        bbox={ 'boxstyle': 'round', 'fc': 'white', 'alpha': 0.9 })
    ax[1].plot((i,i), (0, norm.cdf(i)), color='red', linestyle='--', alpha=0.5)

# annotation on survival zoom
for i in range(6):
    if i < 3:
        text_pos = (i+0.1, 0.5*10**-(i+2)) 
    else:
        text_pos = (i-1.0, 0.5*10**-(i-2))
    fmt='{:.'+f'{i}'+'%}'
    ax[3].annotate(fmt.format(norm.sf(i)),
        xy=(i, norm.sf(i)),
        xytext=text_pos,
        arrowprops={ 'arrowstyle': '->' })
    # vertical lines
    ax[3].plot((i,i), (0, norm.sf(i)), color='red', linestyle='--', alpha=0.5)

plt.show()
