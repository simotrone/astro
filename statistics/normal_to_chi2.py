from scipy.stats import norm, chi2
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches

x = np.linspace(-5, 5, 100)
x_positive = [i for i in x if i >= 0]
# random numbers with gaussian distribution
size = 10_000
r = norm.rvs(size=size)
r2 = r**2

# plots
factor = 1.7
rows,cols=3,2
fig, axes = plt.subplots(rows, cols,
    figsize=(cols*3*factor, rows*2*factor),
    tight_layout=False)
fig.suptitle('Gaussian things II')
ax = axes.flatten()

# row 1
ax[0].plot(x, norm.pdf(x), color='red', ls='-', label='pdf')
ax[0].set_title('Normal pdf + random data')
ax[0].hist(r, density=True, bins=20, histtype='stepfilled', alpha=0.3, label='data')
ax[0].legend(loc='best')

ax[1].plot(x_positive, chi2.pdf(x_positive, df=1), color='darkgreen', ls='-.', label='χ² pdf')
ax[1].set_title('χ² pdf')

# row 2
ax[2].hist(r2, bins=40, histtype='stepfilled', alpha=0.3, label='data²')
ax[2].set_title('Random data squared')

ax[3].plot(x_positive, chi2.pdf(x_positive, df=1), color='darkgreen', ls='-.', label='χ² pdf')
ax[3].set_title('χ² pdf w/ log yaxis')
ax[3].set_yscale('log')

# row 3
ax[4].hist(r2, density=True, bins=40, histtype='stepfilled', alpha=0.3, label='data²')
ax[4].set_title('Data² freq w/ log yaxis')
ax[4].set_yscale('log')
ax[4].set_ylim([1/(size*10), 10])

n, bins, patches = ax[5].hist(r2, density=True, bins=40, histtype='stepfilled', alpha=0.3, label='data²')
x_chi2 = np.linspace(0, int(np.max(bins))+1, 100)
ax[5].set_title('Data² freq and χ² dist')
ax[5].set_yscale('log')
ax[5].set_ylim([1/(size*10), 10])
ax[5].plot(x_chi2, chi2.pdf(x_chi2, df=1), color='darkgreen', ls='-.', label='χ² pdf', alpha=0.7)
ax[5].legend(loc='best')

# arrow patches
# Thanks to: https://www.cilyan.org/blog/2016/01/23/matplotlib-draw-between-subplots/
# Create the arrow
# 1. Get transformation operators for axis and figure
# Transformation tutorial: https://matplotlib.org/tutorials/advanced/transforms_tutorial.html
axtr = []
for i, j in enumerate(ax):
    axtr.append(ax[i].transAxes)   # Axis #i -> Display
figtr = fig.transFigure.inverted() # Display -> Figure

# 2. Transform arrow start point from axis 0 to figure coordinates
# 3. Transform arrow end point from axis 1 to figure coordinates
# The function take (begin_ax_index, end_ax_index, begin_coord, end_coord)
def get_arrow_begin_end(i, j, coord_i=(0.95, 0.1), coord_j=(0.95, 0.9)):
    begin = figtr.transform(axtr[i].transform(coord_i))
    end = figtr.transform(axtr[j].transform(coord_j))
    return (begin, end)

be_02 = get_arrow_begin_end(0, 2)
be_24 = get_arrow_begin_end(2, 4)
be_13 = get_arrow_begin_end(1, 3, coord_i=(0.1, 0.1), coord_j=(0.1, 0.9))
be_35 = get_arrow_begin_end(3, 5, coord_i=(0.1, 0.1), coord_j=(0.1, 0.9))
be_45 = get_arrow_begin_end(4, 5, coord_i=(0.95, 0.85), coord_j=(0.05, 0.85))

# 4. Create the patch
arrows = []
for p in [be_02, be_24, be_45]:
    ptB, ptE = p
    arrows.append(
        matplotlib.patches.FancyArrowPatch(
            ptB, ptE,
            transform=fig.transFigure,                  # Place arrow in figure coord system
            fc = "b", arrowstyle='simple', alpha = 0.3,
            mutation_scale = 40.
        )
    )

for p in [ be_13, be_35 ]:
    ptB, ptE = p
    arrows.append(
        matplotlib.patches.FancyArrowPatch(
            ptB, ptE,
            transform=fig.transFigure,                  # Place arrow in figure coord system
            fc = "darkgreen", arrowstyle='simple', alpha = 0.3,
            mutation_scale = 40.
        )
    )


# 5. Add patch to list of objects to draw onto the figure
for a in arrows:
    fig.patches.append(a)

plt.show()
