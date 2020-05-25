import lib.significance as li_ma

def show_info(s):
    fmt  = 'Simulation {name}\n'
    fmt += '   t_on {t_on:6d}\n'
    fmt += '      α {alpha:12.5f}\n'
    fmt += '   N_on {N_on:8.1f}±{N_on_err:5.1f}\n'
    fmt += '    N_B {N_B:8.1f}±{N_B_err:5.1f}\n'
    fmt += '    N_S {N_S:8.1f}±{N_S_err:5.1f}\n'
    fmt += '  S_(5) {S_(5):8.1f}\n'
    fmt += '  S_(9) {S_(9):8.1f}\n'
    fmt += ' S_(17) {S_(17):8.1f}\n'
    print(fmt.format(**s))

def elaboration(sims, bkg):
    N_off_rel_err = bkg['N_off_err'] / bkg['N_off']
    for s in sims:
        s['alpha'] = s['t_on'] / bkg['t_off']

        s['N_B'] = s['alpha'] * bkg['N_off']
        s['N_B_err'] = s['N_B'] * N_off_rel_err

        s['N_S'] = s['N_on'] - s['N_B']
        s['N_S_err'] = s['N_on_err'] + s['N_B_err']

        # significance equation (5)
        s['S_(5)'] = li_ma.eq_5(s['N_on'], bkg['N_off'], s['alpha'])
        s['S_(9)'] = li_ma.eq_9(s['N_on'], bkg['N_off'], s['alpha'])
        s['S_(17)'] = li_ma.eq_17(s['N_on'], bkg['N_off'], s['alpha'])
    return sims
    
def main():
    # background data
    only_background = { 't_off': 600, 'N_off': 1147, 'N_off_err': 34 }
    # simulations data
    simulations = [
        {'name': '#1', 't_on': 1800, 'N_on': 13814, 'N_on_err': 118},
        {'name': '#2', 't_on':  100, 'N_on':   747, 'N_on_err':  27},
        {'name': '#3', 't_on':   10, 'N_on':    76, 'N_on_err':   9},
    ]

    out = elaboration(simulations, only_background)
    for s in out:
        show_info(s)

main()
