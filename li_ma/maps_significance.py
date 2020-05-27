import argparse
import matplotlib.pyplot as plt
import numpy as np
import lib.significance as li_ma

def significances_matrix(alpha, ons, offs, fn):
    """
    Provide a significances matrix for each on/off values.
    """
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
        # two options to manage the range in colorbar
        if False:
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
    alphas = [1/60, 0.20, 0.5, 1, 2, 10]
    N_on_range  = range(1, 51)
    N_off_range = range(1, 51)
    functions = [
        { 'name': 'eq. 5',  'fn': li_ma.eq_5, },
        { 'name': 'eq. 9',  'fn': li_ma.eq_9, },
        { 'name': 'eq. 17', 'fn': li_ma.eq_17 },
    ]

    n_fns = len(functions)
    n_alphas = len(alphas)

    data = []
    for i, fn in enumerate(functions):
        for j, a in enumerate(alphas):
            data.append({
                'name': fn['name'],
                'significances': significances_matrix(a, N_on_range, N_off_range, fn['fn']),
                'alpha': a,
                'axes': (i, j),
            })

    # plots
    fig, axes = plt.subplots(
        figsize=(n_alphas*3, n_fns*3), # +inches for colorbar
        nrows=n_fns, ncols=n_alphas,
        # gridspec_kw={ 'wspace': 0.05, }
        tight_layout=not OPTS.share_scale,
    )
    fig.suptitle('Significance map plots')
 
    prepare_subplots(
        axes,
        data,
        on_ticks=find_ticks(N_on_range),
        off_ticks=find_ticks(N_off_range),
        colorbar=OPTS.share_scale,
    )
    if OPTS.share_scale:
        fig.colorbar(axes[0][0].get_images()[0], ax=axes, location='top', fraction=0.07)
    else:
        for ax in axes.flatten():
            fig.colorbar(ax.get_images()[0], ax=ax, fraction=0.05)
 
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="ImageMap significance equations from Li&Ma paper")
    parser.add_argument('-s', '--share-scale', action='store_true', help="impose the same significance scale for all the plots.")
    OPTS = parser.parse_args()
    main()

