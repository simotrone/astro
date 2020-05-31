import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from scipy.stats import poisson
# from os import path

def prepare_dist(ax, curves, dist, xlimit=(0,100), title=None):
    x = range(xlimit[1])
    ax.set_xlim(xlimit)
    for c in curves:
        y = dist(x, c['expected'])
        ax.plot(x, y, color=c['color'], linewidth=2, alpha=0.5, label=f'μ: {c["expected"]}')
    if title:
        ax.set_title(title)
    ax.set_xlabel('n')
    ax.set_ylabel('prob')
    ax.legend(loc='best')


def show_static_stats(axes):
    config = [
        { 'expected':  2, 'color': 'blue'  },
        { 'expected': 10, 'color': 'red'   },
        { 'expected': 50, 'color': 'green' },
    ]

    prepare_dist(axes[0], config, poisson.pmf, xlimit=(-5, 80), title='Probability mass function')
    prepare_dist(axes[1], config, poisson.cdf, xlimit=(-1, 80), title='Cumulative distribution function')


def data_gen_setup_simple(mu=10,limit=20):
    def data_gen():
        rv = poisson(mu)
        for i in range(limit):
            yield (i, rv.pmf(i), rv.cdf(i), rv.sf(i))
    return data_gen

# as data_gen_setup_simple, but return arrays of value
# get the [-1] element to have data_gen_setup_simple
def data_setup(mu=10, limit=20):
    rv = poisson(mu)
    step, pmf, cdf, sf = [], [], [], []
    for i in range(limit):
        step.append(i)
        pmf.append(rv.pmf(i))
        cdf.append(rv.cdf(i))
        sf.append(rv.sf(i))
    return (step, pmf, cdf, sf)


def init_setup(axes, artists, xlim=(0,1), ylim=[(0,1),(0,1)],
               xticks=None, xticklabels=None, title=None):
    def init():
        for i, ax in enumerate(axes):
            ax.clear()
            ax.set_xlabel('n')
            ax.set_ylabel('prob')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim[i])
            if title:
                ax.set_title(title[i])
            if xticks:
                ax.set_xticks(xticks)
            if xticklabels:
                ax.set_xticklabels(xticklabels)

        for a in artists:
            if callable(getattr(a, 'set_data', None)):
                a.set_data([],[])
            elif callable(getattr(a, 'set_text', None)):
                a.set_text(f'  n: {"- ":>7s}\npmf: {"- ":>7s}\ncdf: {"- ":>7s}\n sf: {"- ":>7s}')
                a.set_bbox({ 'boxstyle': 'round', 'fc': 'white', 'alpha': 0.9 })
            else:
                raise Exception ('artist', a, 'no data reset?')
        return artists
    return init


# https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.animation.FuncAnimation.html#matplotlib.animation.FuncAnimation
def show_animated_numbers(axes, figure):
    delay_ms = 25
    expected = 14
    n_max = 33

    data = data_setup(expected, n_max)
    main_ax = axes[0]
    support_ax = axes[1]

    # changing artists
    line, = main_ax.plot([], [], color='red', lw=2, alpha=0.7)
    vbar, = main_ax.plot([], [], color='blue', ls=':')
    mu_vline, = main_ax.plot([], [], color="darkgreen", ls='--', label='expected')
    tbox = main_ax.text(0.95, 0.95, s='',
        transform=main_ax.transAxes, ha='right', va='top',
        fontfamily='monospace')
    patch = patches.Rectangle((0,0), 0, 0.01, fc='darkgreen', alpha=0.3)
    main_ax.add_patch(patch)
    sfbox = main_ax.text(n_max-1, -0.010, s='-', fontfamily='monospace', ha='right', va='center')

    cline, = support_ax.plot([], [], color='blue', lw=2, alpha=0.7)
    cvbar, = support_ax.plot([], [], color='blue', lw=2, alpha=0.3)
    mu_cvline, = support_ax.plot([], [], color='darkgreen', ls='--', label='expected')
    ctbox = main_ax.text(0.05, 0.95, s='',
        transform=support_ax.transAxes, ha='left', va='top',
        fontfamily='monospace')

    # plot initialization
    init = init_setup(
        [main_ax, support_ax],
        artists=[line, vbar, mu_vline, tbox, cline, cvbar, mu_cvline, ctbox],
        xlim=(-1, n_max),
        ylim=[(-0.021, max(data[1])+0.05), (-0.17, 1.3)],
        xticks=range(n_max),
        xticklabels=[ i if not i % 5 else '' for i in range(n_max)],
        title=[
            f'Probability Mass Function (μ: {expected})',
            f'Cumulative Distribution F. (μ: {expected})' 
        ]
    )

    def animate(i, n_max, data): #step, pmf, cdf, sf):
        n = i % n_max
        step, pmf, cdf, sf = data

        # always in loop
        vbar.set_data([n,n], [0, pmf[n]])

        # only in the first loop
        if i < n_max:
            line.set_data(step[0:n+1], pmf[0:n+1])
        if i < n_max and n == expected:
            mu_vline.set_data([n,n], [0, pmf[n]])

        # text box changes info through the loops
        tbox_fmt  = f'n: {n:6d} \n'
        tbox_fmt += f'pmf: {pmf[n]:7.2%}\n'
        if i < n_max:
            tbox_fmt += f'cdf: {"- ":>7s}\n'
            tbox_fmt += f'sf: {"- ":>7s}'
        else:
            tbox_fmt += f'cdf: {cdf[n]:7.2%}\n'
            tbox_fmt += f'sf: {sf[n]:7.2%}'
        tbox.set_text(tbox_fmt)

        # only in the second loop
        if i >= n_max:
            zeros = [0] * (n+1)
            main_ax.collections.clear()
            main_ax.fill_between(step[0:n+1], zeros, pmf[0:n+1], color='blue', alpha=0.3)
            main_ax.figure.canvas.draw()

        if i >= n_max:
            patch.set_width(n_max-1-n)
            patch.set_xy([n, -0.015])
            sfbox.set_text(f'survival:{sf[n]:7.2%} ')
        else:
            sfbox.set_text('')

        if i >= n_max:
            cline.set_data(step[0:n+1], cdf[0:n+1])
            cvbar.set_data([n,n], [0, cdf[n]])
        if i >= n_max and n == expected:
            mu_cvline.set_data([n,n], [0, cdf[n]])

        # cumulative textbox
        if i >= n_max:
            ctbox_fmt  = f'  n: {n:6d} \n'
            ctbox_fmt += f'pmf: {"- ":>7s}\n'
            ctbox_fmt += f'cdf: {cdf[n]:7.2%}\n'
            ctbox_fmt += f' sf: {sf[n]:7.2%}'
            ctbox.set_text(ctbox_fmt)

        return line, vbar, mu_vline, tbox, cline, cvbar, mu_cvline, ctbox, patch, sfbox

    anim = animation.FuncAnimation(figure, animate, fargs=(n_max, data),
        init_func=init, frames=n_max*2,
        interval=delay_ms, repeat_delay=500, blit=True)

    return anim



def main():
    factor=1.0
    rows=2
    cols=2

    fig, axes = plt.subplots(rows, cols,
        figsize=(cols*6.4*factor, rows*4.8*factor))
    fig.suptitle('Poisson stats')

    show_static_stats(axes[0])
    anim = show_animated_numbers(axes[1], fig)

    # TODO: animation save
    # filename = 'poisson_plots.mp4'
    # if not path.exists(filename):
    #     ...
    # anim.save('poisson_plots.mp4', writer='ffmpeg', codec='libx264', fps=15)
    # anim.save('poisson_plots.gif', writer='imagemagick', fps=15)

    plt.show()


if __name__ == '__main__':
    main()


