import zlib # need to avoid a "core seg fault" when using matplotlib
import cscripts
import ctools
import gammalib
import json
import math
import matplotlib.pyplot as plt
import os
import statistics
import sys
from scipy.stats import norm
# version 2.1
# Example:
#   MAX_SEED=10 python exec.py
#   ds9 cntcube.fits -region onoff_on.reg 
# TODO
#  . see flux eval with gammapy lib

OBJ = { "ra": 83.6331, "dec": 22.0145, "model": "crab.xml", "name": "Crab" }
ENERGY = { "min": 0.03, "max": 150.0 }
TIME_SLOTS = [1800, 600, 100, 60, 30, 10, 5]
SEED_RANGE = [1]
BINNING = False # output binning file take a lot of MB
POINTINGS = [
    {
        "ra": OBJ["ra"],
        "dec": OBJ["dec"]+0.5,
        "name": "dec_0.5",
    },
    {
        "ra": OBJ["ra"],
        "dec": OBJ["dec"]+1.0,
        "name": "dec_1.0",
    },
    {
        "ra": OBJ["ra"]+0.5,
        "dec": OBJ["dec"],
        "name": "ra_0.5",
    },
    {
        "ra": OBJ["ra"]+1.0,
        "dec": OBJ["dec"]+1.0,
        "name": "ra_1.0_dec_1.0",
    },
    {
        "ra": OBJ["ra"]+2.0,
        "dec": OBJ["dec"]+2.0,
        "name": "ra_2.0_dec_2.0",
    },
]

def create_dir(d):
    ret = False
    if not os.path.isdir(d):
        os.mkdir(d)
        ret = True
    return ret
    
def simulation_run(output_events_file, ra, dec, tmax=1800, seed=1, force=0):
    log_file = os.path.join(os.path.dirname(output_events_file), "ctobssim.log")
    if not os.path.isfile(output_events_file) or force == 1:
        sim = ctools.ctobssim()
        sim.clear()
        sim["seed"]  = seed
        sim["ra"]    = ra
        sim["dec"]   = dec
        sim["rad"]   = 5.0
        sim["tmin"]  = 0
        sim["tmax"]  = tmax
        sim["emin"]  = ENERGY["min"]
        sim["emax"]  = ENERGY["max"]
        sim["caldb"] = "prod2"
        sim["irf"]   = "South_0.5h"
        sim["inmodel"]   = OBJ["model"]
        sim["outevents"] = output_events_file
        sim["logfile"]   = log_file
        sim.logFileOpen()
        sim.run()
        sim.save()
        sys.stderr.write("File '%s' created.\n" % output_events_file)
    return True

# taglia gli eventi solo su tmax
def select_events(events_file, selected_file, tmax=100, force=0):
    if tmax < 1:
        print("need tmax > 0 to select events")
        exit(1)
    log_file = os.path.join(os.path.dirname(selected_file), "ctselect_"+str(tmax)+".log")
    if not os.path.isfile(selected_file) or force == 1:
        select = ctools.ctselect()
        select.clear()
        select["inobs"] = events_file
        select["rad"] = "INDEF"
        select["tmin"] = 0
        select["tmax"] = tmax
        select["emin"] = "INDEF"
        select["emax"] = "INDEF"
        select["outobs"] = selected_file
        select["logfile"] = log_file
        select.logFileOpen()
        select.run()
        select.save()
        sys.stderr.write("File '%s' created.\n" % selected_file)
    return True

def ctbinning_run(ev_file, force=0):
    working_dir = os.path.dirname(ev_file)
    cnt_file = os.path.join(working_dir, "cntcube.fits")
    log_file = os.path.join(working_dir, "ctbin.log")
    if not os.path.isfile(cnt_file) or force == 1:
        binning = ctools.ctbin()
        binning.clear()
        binning["inobs"]   = ev_file
        binning["outobs"]  = cnt_file
        binning["xref"]    = OBJ["ra"]
        binning["yref"]    = OBJ["dec"]
        binning["ebinalg"] = "LOG"
        binning["emin"]    = ENERGY["min"]
        binning["emax"]    = ENERGY["max"]
        binning["enumbins"] = 20
        binning["nxpix"]    = 500
        binning["nypix"]    = 500
        binning["binsz"]    = 0.01
        binning["coordsys"] = "CEL"
        binning["proj"]     = "CAR"
        binning["logfile"]  = log_file
        binning.logFileOpen()
        binning.run()
        binning.save()
        sys.stderr.write("File '%s' created.\n" % cnt_file)
    return cnt_file

# http://cta.irap.omp.eu/ctools/users/user_manual/classical.html
def csphagen_run(ev_file, force=0):
    working_dir = os.path.dirname(ev_file)
    log_file = os.path.join(working_dir, "csphagen.log")
    output_model_file = os.path.join(working_dir, "onoff_result.xml")
    output_obs_def = os.path.join(working_dir, "onoff_obs.xml")
    if not os.path.isfile(output_obs_def) or not os.path.isfile(output_model_file) or force == 1:
        phagen = cscripts.csphagen()
        phagen.clear()
        phagen["inobs"] = ev_file
        phagen["inmodel"] = OBJ["model"]
        phagen["srcname"] = OBJ["name"]
        phagen["caldb"] = "prod2"
        phagen["irf"]   = "South_0.5h"
        phagen["ebinalg"]  = "LIN"
        phagen["emin"]     = ENERGY["min"]
        phagen["emax"]     = ENERGY["max"]
        phagen["enumbins"] = 1
        phagen["coordsys"] = "CEL"
        phagen["ra"]    = OBJ["ra"]
        phagen["dec"]   = OBJ["dec"]
        phagen["rad"]   = 0.2
        phagen["stack"] = False
        phagen["bkgmethod"] = "REFLECTED"
        phagen["outobs"]   = output_obs_def
        phagen["outmodel"] = output_model_file
        phagen["prefix"]   = os.path.join(working_dir, "onoff")
        phagen["logfile"]  = log_file
        phagen.logFileOpen()
        phagen.run()
        phagen.save()
        sys.stderr.write("File '%s' created.\n" % output_obs_def)
        sys.stderr.write("File '%s' created.\n" % output_model_file)
    return { "obs": output_obs_def, "model": output_model_file }

def ctlike_run(observation_file, input_model, force=0):
    working_dir = os.path.dirname(observation_file)
    result_file = os.path.join(working_dir, "ml_result.xml")
    log_file    = os.path.join(working_dir, "ctlike.log")
    if not os.path.isfile(result_file) or force == 1:
        like = ctools.ctlike()
        like.clear()
        like["inobs"]    = observation_file
        like["inmodel"]  = input_model
        like["outmodel"] = result_file
        like["logfile"]  = log_file
        like.logFileOpen()
        like.run()
        like.save()
        sys.stderr.write("File '%s' created.\n" % result_file)
    return result_file

def compute_fluxes(spectral_model, emin_tev, emax_tev):
    e_min = gammalib.GEnergy(emin_tev, "TeV")
    e_max = gammalib.GEnergy(emax_tev, "TeV")
    ret = {
        "flux":  spectral_model.flux(e_min, e_max),  # ph/cm²/s
        "eflux": spectral_model.eflux(e_min, e_max), # erg/cm²/s
    }
    return ret

def get_model_info(model_file):
    models = gammalib.GModels(model_file)
    spectral_model = models["Crab"].spectral()
    fluxes = compute_fluxes(spectral_model, ENERGY["min"], ENERGY["max"])
    return {
        "model_file": model_file,
        "working_dir": os.path.dirname(model_file),
        "index": {
            "value": models["Crab"]["Index"].value(),
            "error": models["Crab"]["Index"].error(),
        },
        "prefactor": {
            "value": models["Crab"]["Prefactor"].value(),
            "error": models["Crab"]["Prefactor"].error(),
        },
        "pivot": {
            "value": models["Crab"]["PivotEnergy"].value(),
            "error": models["Crab"]["PivotEnergy"].error(),
        },
        "ts": models["Crab"].ts(),
        "flux": fluxes["flux"], # ph/cm²/s
        "eflux": fluxes["eflux"], # erg/cm²/s
    }

# verified with https://docs.gammapy.org/0.8/stats/index.html#li-ma-significance
def li_ma (n_on, n_off, alpha):
    if n_on <= 0 or n_off <= 0 or alpha == 0:
        return None
    fc = (1 + alpha) / alpha
    fb = n_on / (n_on + n_off)
    f  = fc * fb
    gc = 1 + alpha
    gb = n_off / (n_on + n_off)
    g  = gc * gb
    first  = n_on * math.log(f)
    second = n_off * math.log(g)
    fullb   = first + second
    return math.sqrt(2) * math.sqrt(fullb)

def get_onoff_info(on_off_obs_file):
    obs = gammalib.GObservations(on_off_obs_file)
    onoff_obs = obs[0] # GCTAOnOffObservation
    # Pulse Height Analyzer class.
    # This class implements a Pulse Height Analyzer (PHA) spectrum that is used
    # as data container for an XSPEC analysis. A PHA spectrum is a vector that
    # provides the number of measured counts as function of the channel number.
    pha_on  = onoff_obs.on_spec()
    pha_off = onoff_obs.off_spec()
    spectrum_bins = pha_on.size()
    if spectrum_bins != 1:
        print("spectrum bins are more then expected. Need to change something in code to manage them")
    on_count  = pha_on.counts()
    off_count = pha_off.counts()
    alpha_val = pha_on.backscal(spectrum_bins-1) # index 0 = spectrum_bins -1 
    li_ma_val = li_ma(on_count, off_count, alpha_val)
    return {
        "N_on": on_count,
        "N_off": off_count,
        "N_s": (on_count - alpha_val * off_count),
        "alpha": alpha_val,
        "li_ma": li_ma_val,
    }

# db_ref[0] example:
# {'dec': 24.0145,
#  'name': 'ra_2.0_dec_2.0',
#  'event_file': 'ra_2.0_dec_2.0/selected_tmax_5/events_tmax_5.fits',
#  'model_file': 'ra_2.0_dec_2.0/selected_tmax_5/ml_result.xml',
#  'onoff_obs': 'ra_2.0_dec_2.0/selected_tmax_5/onoff_obs.xml',
#  'ra': 85.6331,
#  'seed': 1,
#  'tmax': 5}]
def data_summary(all_data_info):
    o = {}
    for i in all_data_info:
        key = i["name"] +"_"+ str(i["tmax"])
        if key not in o:
            o[key] = {
                "name": i["name"],
                "tmax": i["tmax"],
                "ra":   i["ra"],
                "dec":  i["dec"],
                "event_file": i["event_file"],
                "model_file": i["model_file"],
                "onoff_obs":  i["onoff_obs"],
                "seed_array":  [],
                "flux_array":  [],
                "eflux_array": [],
                "ts_array":      [],
                # "sqrt_ts_array": [],
                "li_ma_array": [],
                "N_on_array":  [],
                "N_off_array": [],
                "N_s_array":   [],
                "alpha_array": [],
            }
        model = get_model_info(i["model_file"])
        onoff = get_onoff_info(i["onoff_obs"])

        o[key]["seed_array"].append(i["seed"])
        o[key]["flux_array"].append(model["flux"])
        o[key]["eflux_array"].append(model["eflux"])
        o[key]["ts_array"].append(model["ts"]) # need ts to 0?
        # o[key]["sqrt_ts_array"].append(math.sqrt(model["ts"]) if model["ts"] >= 0 else 0)
        o[key]["li_ma_array"].append(onoff["li_ma"] if onoff["li_ma"] is not None else 0)
        o[key]["N_on_array"].append(onoff["N_on"])
        o[key]["N_off_array"].append(onoff["N_off"])
        o[key]["N_s_array"].append(onoff["N_s"])
        o[key]["alpha_array"].append(onoff["alpha"]) # useless, just for check

        if model["ts"] < 0:
            sys.stderr.write("%15s (%d on, %2d off, %3d seed, %4d tmax): Negative ts %.2f\n" % (i["name"], onoff["N_on"], onoff["N_off"], i["seed"], i["tmax"], model["ts"]))
        elif onoff["li_ma"] is None:
            sys.stderr.write("%15s (%d on, %2d off, %3d seed, %4d tmax): Cannot calculate Li&Ma\n" % (i["name"], onoff["N_on"], onoff["N_off"], i["seed"], i["tmax"]))

    return o

def data_summary_is_ok(data):
    # check size
    if len(data) != len(POINTINGS)*len(TIME_SLOTS):
        sys.stderr.write("Data summary length is %d and should be %d (pointings x time_slots)\n" % (len(data), len(POINTINGS)*len(TIME_SLOTS)))
        return False
    for k in data:
        # check array data
        for sk in data[k]:
            if type(data[k][sk]) != type([]):
                continue
            if len(data[k][sk]) == len(SEED_RANGE):
                continue
            sys.stderr.write("not enough data for '%s'\n" % k)
            sys.stderr.write("  key '%s' has %d values and should be %d\n" % (sk, len(data[k][sk]), len(SEED_RANGE)))
            return False
    return True

def array_stats(arr):
    stat = {
        "n": len(arr),
        "mean": statistics.mean(arr),
        "stdev": statistics.pstdev(arr),
        "median": statistics.median(arr),
    }
    # same numbers of: mean, stdev = norm.fit(arr)
    return stat

def print_data_summary(data):
    fields = [
      #   h_format,      v_format,             title,             sub_t
        [   "%15s",        "%15s",          "fs ref",      "==========", ],
        [   "%10s",      "%10.4f",              "RA",              "==", ],
        [   "%10s",      "%10.4f",             "Dec",             "===", ],
        [    "%6s",         "%6d",            "tmax",            "====", ],
        [    "%6s",         "%6d",           "seeds",           "=====", ],
        [   "%18s", "%9.2e±%8.2e", "flux [ph/cm²/s]", "===============", ],
        [   "%16s", "%9.2f±%6.2f",              "TS",              "==", ],
      # [    "%6s",       "%6.2f",             "√TS",             "===", ],
        [   "%15s", "%8.2f±%6.2f",            "N_on",            "====", ],
        [   "%15s", "%8.2f±%6.2f",           "N_off",           "=====", ],
        [   "%15s", "%8.2f±%6.2f",             "N_s",             "===", ],
        [   "%11s", "%6.2f±%4.2f",           "Li&Ma",           "=====", ],
        [    "%7s",       "%7.4f",           "alpha",           "=====", ],
    ]

    header_fmt = " ".join([r[0] for r in fields]) # headers format
    values_fmt = " ".join([r[1] for r in fields]) # values format
    print(header_fmt % tuple([r[2] for r in fields])) # titles
    print(header_fmt % tuple([r[3] for r in fields])) # sub_titles separator
    for d in sorted(data, key=lambda i: (-1*i["tmax"], i["ra"], i["dec"])):
        n_seeds = len(d["seed_array"])
        flux_m = array_stats(d["flux_array"])
        ts_m   = array_stats(d["ts_array"])
        # sqrt_ts_m = array_stats(d["sqrt_ts_array"])
        N_on_m  = array_stats(d["N_on_array"])
        N_off_m = array_stats(d["N_off_array"])
        N_s_m   = array_stats(d["N_s_array"])
        li_ma_m = array_stats(d["li_ma_array"])
        alpha_m = array_stats(d["alpha_array"]) # useless
        print(values_fmt % (d["name"], d["ra"], d["dec"], d["tmax"], n_seeds,
            flux_m["mean"], flux_m["stdev"],
            ts_m["mean"], ts_m["stdev"],
            # sqrt_ts_m["mean"],
            N_on_m["mean"],  N_on_m["stdev"],
            N_off_m["mean"], N_off_m["stdev"],
            N_s_m["mean"],   N_s_m["stdev"],
            li_ma_m["mean"], li_ma_m["stdev"],
            alpha_m["mean"]))

def dynamic_bin_number(arr, max_val=None):
    n = max(arr)-min(arr)
    if n < 1:
        n = len(arr)
    if max_val is not None and n > max_val:
        n = max_val
    return int(n)

# seaborn graph with distplot. Same data, same gaussian loc/scale
#   import seaborn as sns, numpy as np
#   print(array_stats(d["ts_array"]))
#   print(norm.fit(d["ts_array"]))
#   sns.distplot(d["ts_array"], bins=50, kde=True, fit=norm, norm_hist=False)# , norm_hist=True) #array, bins=n_bins, fit=norm, norm_hist=True
def create_hist(ax, data_arr, xlabel=None, color="blue", dyn_bins=False, exp_fmt=False, density_flag=True):
    n_bins = 50
    if dyn_bins:
        n_bins = dynamic_bin_number(data_arr)
    stats = array_stats(data_arr)

    if density_flag:
        counts, bins, patches = ax.hist(data_arr, bins=n_bins, alpha=0.5, edgecolor=color, color=color, density=True)
        ax.plot(bins, norm.pdf(bins, stats["mean"], stats["stdev"]), color="grey", linestyle="--", alpha=0.9)
    else:
        counts, bins, patches = ax.hist(data_arr, bins=n_bins, alpha=0.5, edgecolor=color, color=color)

    ax.axvline(stats["mean"],   color="grey", linestyle="--", alpha=0.9)
    ax.axvline(stats["median"], color="black",  linestyle="-.", alpha=0.5)

    legend_fmt = 'bins=%d\nmean=%.2f\nmedian=%.2f\nσ=%.2f'
    if exp_fmt:
        legend_fmt = 'bins=%d\nmean=%.2e\nmedian=%.2e\nσ=%.2e'
    ax.text(x=0.65, y=0.70, s=legend_fmt % (n_bins, stats["mean"], stats["median"], stats["stdev"]), transform=ax.transAxes)

    # if int is True:
    #     from matplotlib.ticker import FormatStrFormatter
    #     bin_centers = [ (bins[i]+bins[i+1])*0.5 for i in range(len(bins)-1) ]
    #     ax.set_xticks(bin_centers)
    #     #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    #     ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    if xlabel is not None:
        ax.set_xlabel(xlabel)

def plot_data_summary(data):
    rows_num=3
    for d in data:
        # figsize in inch (width, height)
        fig, ax = plt.subplots(nrows=rows_num, ncols=2, figsize=(9, rows_num*3))
        fig.suptitle(d["name"]+" "+str(d["tmax"]), va="top", ha="center")

        create_hist(ax[0][0], d["ts_array"],    color="blue",   xlabel="TS")
        create_hist(ax[1][0], d["li_ma_array"], color="cyan",   xlabel="Li&Ma")
        create_hist(ax[2][0], d["flux_array"],  color="green",  xlabel="Flux [ph/cm²/s]", exp_fmt=True)
        create_hist(ax[0][1], d["N_on_array"],  color="orange", xlabel="N on",     dyn_bins=True)
        create_hist(ax[1][1], d["N_off_array"], color="red",    xlabel="N off",    dyn_bins=True)
        create_hist(ax[2][1], d["N_s_array"],   color="yellow", xlabel="N signal", dyn_bins=True)

        # Important: first tight_layout(), after adjust for the title
        fig.tight_layout()
        fig.subplots_adjust(top=0.95)
        # plt.show()
        plt.savefig("%s_%04d.png" % (d["name"], d["tmax"]), format="png")
        plt.close()

    return None

# from Jurgen's ctools example
# Ex: plot_spectrum("crab.xml")
def plot_spectrum(model_file, emin=0.03, emax=150.0, enumbins=100):
    models = gammalib.GModels(model_file)
    spectral_model = models["Crab"].spectral()
    e_min = gammalib.GEnergy(emin, "TeV")
    e_max = gammalib.GEnergy(emax, "TeV")
    ebounds = gammalib.GEbounds(enumbins, e_min, e_max)

    x = []
    y = []
    for i in range(enumbins):
        energy = ebounds.elogmean(i)
        value = spectral_model.eval(energy)
        x.append(energy.TeV())
        y.append(value)

    plt.figure()
    plt.title(models["Crab"].name()+" ("+spectral_model.type()+")")
    plt.grid()
    plt.loglog(x, y, color="red")
    plt.xlabel("Energy (TeV)")
    plt.ylabel("dN/dE ph/cm²/s/MeV")
    plt.show()

# data for every time_slot in output reference
# { event_file: , tmax: , ra: , dec: , name: , seed: }
def create_data(pnt, time_slots=[1800], output=[], seed=1, force=0):
    # manage the first time slot
    starting_tmax = time_slots[0]
    if create_dir(pnt["name"]):
        sys.stderr.write("Dir '%s' created.\n" % pnt["name"])
    working_dir = os.path.join(pnt["name"], str(seed))
    create_dir(working_dir)
    events_file = os.path.join(working_dir, "events.fits")

    simulation_run(events_file, ra=pnt["ra"], dec=pnt["dec"], tmax=starting_tmax, seed=seed, force=force)
    output.append({ "event_file": events_file, "tmax": starting_tmax, "ra": pnt["ra"], "dec": pnt["dec"], "name": pnt["name"], "seed": seed })

    # manage other timings in time_slots
    for t in time_slots[1:]:
        sel_working_dir = os.path.join(working_dir, "sel_"+str(t))
        create_dir(sel_working_dir)
        sel_out_file = os.path.join(sel_working_dir, "events_"+str(t)+".fits")
        select_events(events_file, sel_out_file, tmax=t, force=force)
        output.append({ "event_file": sel_out_file, "tmax": t, "ra": pnt["ra"], "dec": pnt["dec"], "name": pnt["name"], "seed": seed })

def run(seed):
    db_ref = []
    # create data for each pointings and put info in db_ref
    for pnt in POINTINGS:
        create_data(pnt, time_slots=TIME_SLOTS, output=db_ref, seed=seed, force=0)

    # run csphagen and ctlike for every generated data
    for d in db_ref:
        # ra, dec, tmax, file
        if BINNING: 
            ctbinning_run(d["event_file"], force=0)
        onoff = csphagen_run(d["event_file"], force=0)
        d["onoff_obs"] = onoff["obs"]
        model_file = ctlike_run(observation_file=onoff["obs"], input_model=onoff["model"], force=0)
        d["model_file"] = model_file

    return db_ref

if __name__ == "__main__":
    if "SEED" in os.environ:
        SEED_RANGE=[int(os.environ["SEED"])]

    if "MAX_SEED" in os.environ:
        SEED_RANGE=range(int(os.environ["MAX_SEED"]))

    sys.stderr.write("Seed range: %s\n" % SEED_RANGE)

    dump_fn = "dump_"+str(max(SEED_RANGE)+1)+".json"
    if "DUMP_FILE" in os.environ:
        dump_fn = os.environ["DUMP_FILE"]

    ds = []
    if not os.path.isfile(dump_fn):
        # create a data struct for each seed
        data = []
        for seed in SEED_RANGE:
            r = run(seed=seed)
            data += r
        sys.stderr.write("Data length: %d\n" % len(data))

        ds = data_summary(data)
        if not data_summary_is_ok(ds):
            exit(1)

        sys.stderr.write("Save data to: %s\n" % dump_fn)
        with open(dump_fn, "w") as f:
            json.dump(ds, f)
    else:
        sys.stderr.write("Load data from: %s\n" % dump_fn)
        with open(dump_fn, "r") as f:
            ds = json.load(f)
        if not data_summary_is_ok(ds):
            exit(1)

    print_data_summary(ds.values())
    plot_data_summary(list(ds.values()))
    exit(0)
