import cscripts as cs
import ctools as ct
import gammalib as gl
import math
import sys
from lib.utils import li_ma
import argparse
# PYTHONPATH=path/to/lib python delta_significance.py onoff_obs_list.xml ml_result.xml 

def inspect_onoff_observations(onoff_obs_file):
    oo_obs_list = gl.GObservations(onoff_obs_file)
    # print("Observations:", len(oo_obs_list)) # len = 1
    onoff_obs = oo_obs_list[0]
    # print(onoff_obs)

    on_spectrum  = onoff_obs.on_spec()
    off_spectrum = onoff_obs.off_spec()
    print("ON\n",  on_spectrum, "\nSize: ", on_spectrum.size())
    print("OFF\n", off_spectrum)
    alpha = on_spectrum.backscal(on_spectrum.size()-1)
    on_counts  = on_spectrum.counts()
    off_counts = off_spectrum.counts()
    excess_counts = on_counts - off_counts * alpha
    sign = li_ma(on_counts, off_counts, alpha)
    print('''pha summary:
        on: {}
       off: {}
     alpha: {}
    excess: {}
     li_ma: {}'''.format(on_counts, off_counts, alpha, excess_counts, sign))
    return sign

def inspect_likelihood_model(model):
    ml = gl.GModels(model)
    # print(ml)
    ts = ml[0].ts()
    sign = math.sqrt(ts)
    print('''likelihood:
        TS: {}
       √TS: {}'''.format(ts, sign))
    return sign

if __name__ == '__main__':
    #onoff_obs_file = 'onoff_obs_list.xml'
    parser = argparse.ArgumentParser(description="Compare onoff vs model ts")
    parser.add_argument("onoff_observation_file", help="the xml file with onoff data (pha)")
    parser.add_argument("-m", "--model", help="ctlike output model")
    args = parser.parse_args()
    lima_sign = inspect_onoff_observations(args.onoff_observation_file)
    if args.model:
        like_sign = inspect_likelihood_model(args.model)
        print('Δ significance: {0:.6f}'.format(lima_sign-like_sign))
