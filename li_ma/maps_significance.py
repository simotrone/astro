import math
import matplotlib.pyplot as plt
from matplotlib import colors
# import mpl_toolkits.axisartist.axislines as axislines
import numpy as np
from lib.utils import li_ma

def significance_eq_5(n_on, n_off, alpha):
    num = n_on - alpha * n_off
    den = math.sqrt(n_on + alpha**2 * n_off)
    return num / den

def map_significances(alpha, ons, offs, fn):
    s = []
    for on in ons:
        row = []
        for off in offs:
            row.append(fn(on, off, alpha))
        s.append(row)
    return s

def prepare_subplots(axes, data, on_ticks=None, off_ticks=None, colorbar=False):
    vmin, vmax = None, None
    if colorbar:
        values_buffer = []
        for d in data:
            values_buffer.append(np.min(d['significances']))
            values_buffer.append(np.max(d['significances']))
        if True:
            # data min and max
            vmin, vmax = np.min(values_buffer), np.max(values_buffer)
        else:
            # absolute data limits scale
            val_max = np.max(np.abs(values_buffer))
            vmin, vmax = -1*val_max, val_max

    for d in data:
        i = d['axes'][0] # row index
        j = d['axes'][1] # col index
        ax = axes[i][j]
        ax.imshow(
            d['significances'],
            origin='lower',
            vmin=vmin, vmax=vmax,
            cmap='RdBu',
            #interpolation='none',
        )

        if off_ticks is not None:
            ax.set_xticks(off_ticks['indexes'])
            ax.set_xticklabels(off_ticks['labels'])
        if on_ticks is not None:
            ax.set_yticks(on_ticks['indexes'])
            ax.set_yticklabels(on_ticks['labels'])

        if ax.is_last_row():
            ax.set(xlabel='N off\n'+r'$\bf{alpha: '+f'{d["alpha"]:.3f}'+r'}$')
        else:
            ax.set_xticklabels([])

        if ax.is_first_col():
            ax.set(ylabel=r'$\bf{'+d["name"]+r'}$'+'\nN on')
        else:
            ax.set_yticklabels([])

def find_ticks(values, how_many=6):
    samples = np.linspace(min(values), max(values), how_many, endpoint=True, dtype=int)
    indexes = [i for i, v in enumerate(values) if v in samples]
    return { 'labels': samples, 'indexes': indexes }

def main():
    alphas = [0.0166666, 0.20, 0.25, 1]
    N_ons = range(1, 21)
    N_offs = range(1, 21)
    functions = [
        { 'name': 'eq. 5',   'fn': significance_eq_5, },
        { 'name': 'eq. 17*', 'fn': li_ma,             },
    ]

    n_fns = len(functions)
    n_alphas = len(alphas)

    data = []
    for i, fn in enumerate(functions):
        for j, a in enumerate(alphas):
            data.append({
                'name': fn['name'],
                'significances': map_significances(a, N_ons, N_offs, fn['fn']),
                'alpha': a,
                'axes': (i, j),
            })

    # plots
    fig, axes = plt.subplots(
        figsize=(n_alphas*3, n_fns*3+1), # +inches for colorbar
        nrows=n_fns, ncols=n_alphas,
        # gridspec_kw={ 'wspace': 0.05, }
    )
    fig.suptitle('Significances plots')
 
    colorbar_flag = False
    prepare_subplots(
        axes,
        data,
        on_ticks=find_ticks(N_ons),
        off_ticks=find_ticks(N_offs),
        colorbar=colorbar_flag
    )
    if colorbar_flag:
        fig.colorbar(axes[0][0].get_images()[0], ax=axes, location='top', fraction=0.07)
 
    # fig.tight_layout()
    plt.show()

main()

