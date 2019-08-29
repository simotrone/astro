import cscripts
import ctools
import gammalib
import os

class CToolsWrapper:
    def __init__(self, args, verbosity=0):
        self.verbosity = verbosity
        fields = [
            ['name', 'MANDATORY'],
            ['ra',   'MANDATORY'], # source ra
            ['dec',  'MANDATORY'], # source dec
            ['seed',       1],
            ['rad',        5.0],
            ['energy_min', 0.03],
            ['energy_max', 150.0],
            ['caldb',    'prod2'],
            ['irf', 'South_0.5h'],
        ]
        for f in fields:
            field, default = f
            if default == 'MANDATORY':
                self.__dict__[field] = args[field]
            else:
                self.__dict__[field] = args[field] if field in args and args[field] else default

    def simulation_run(self, model_file, events_file, ra=None, dec=None, time=[0, 1800], energy=[None, None], log_file='ctobssim.log', force=False, save=False):
        sim = ctools.ctobssim()
        sim["inmodel"]   = model_file
        sim["outevents"] = events_file
        sim["seed"]  = self.seed
        sim["ra"]    = ra  if ra  else self.ra
        sim["dec"]   = dec if dec else self.dec
        sim["rad"]   = self.rad
        sim["tmin"]  = float(time[0])
        sim["tmax"]  = float(time[1])
        sim["emin"]  = energy[0] if energy[0] else self.energy_min
        sim["emax"]  = energy[1] if energy[1] else self.energy_max
        sim["caldb"] = self.caldb
        sim["irf"]   = self.irf
        sim["logfile"] = log_file
        if force or not os.path.isfile(events_file):
            sim.run()
            sim.logFileOpen()
        else:
            container = gammalib.GObservations()
            gcta_obs = gammalib.GCTAObservation(events_file)
            container.append(gcta_obs)
            sim.obs(container)
        saved = False
        if (save and force) or (save and not os.path.isfile(events_file)):
            sim.save()
            saved = True
        if saved and self.verbosity > 1:
            print("Events file '{}' saved. time [{}-{}]".format(sim["outevents"].value(), tstart, tstop), file=sys.stderr)
        return sim

    def csphagen_run(self, input_obs_list, input_model, source_rad=0.2, output_obs_list='onoff_obs.xml', output_model='onoff_result.xml', log_file='csphagen.log', prefix='onoff', force=False, save=False):
        phagen = cscripts.csphagen()
        if isinstance(input_obs_list, gammalib.GObservations):
            phagen.obs(input_obs_list)
        elif os.path.isfile(input_obs_list):
            # observations list from file
            phagen["inobs"] = input_obs_list
        else:
            raise Exception('Cannot understand input obs list for csphagen')
        phagen["inmodel"] = input_model
        phagen["srcname"] = self.name
        phagen["caldb"]   = self.caldb
        phagen["irf"]     = self.irf
        phagen["ebinalg"]  = "LIN"
        phagen["emin"]     = self.energy_min
        phagen["emax"]     = self.energy_max
        phagen["enumbins"] = 1
        phagen["coordsys"] = "CEL"
        phagen["ra"]    = self.ra
        phagen["dec"]   = self.dec
        phagen["rad"]   = source_rad
        phagen["stack"] = True
        phagen["bkgmethod"] = "REFLECTED"
        phagen["outobs"]   = output_obs_list
        phagen["outmodel"] = output_model
        phagen["prefix"]   = prefix
        phagen["logfile"]  = log_file
        if force or not os.path.isfile(output_obs_list) or not os.path.isfile(output_model):
            phagen.run()
            phagen.logFileOpen()
        elif os.path.isfile(output_obs_list) and os.path.isfile(output_model):
            onoff_obs = gammalib.GObservations(output_obs_list)
            onoff_obs.models(gammalib.GModels(output_model))
            phagen = cscripts.csphagen(onoff_obs)
        else:
            raise Exception("Cannot proceed with csphagen")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_obs_list)) or (save and not os.path.isfile(output_model)):
            phagen.save()
            saved = True
        if saved and self.verbosity > 1:
            print("Files {}, {} created.".format(output_obs_list, output_model), file=sys.stderr)
        return phagen

    def ctlike_run(self, input_obs_list, input_models=None, output_models='ml_result.xml', log_file='ctlike.log', force=False, save=False):
        like = ctools.ctlike()
        if isinstance(input_obs_list, gammalib.GObservations):
            like.obs(input_obs_list)
        elif os.path.isfile(input_obs_list) and os.path.isfile(input_models):
            # observations list from file
            like["inobs"] = input_obs_list
            like["inmodel"] = input_models
        else:
            raise Exception('Cannot understand input obs list for ctlike')
        like["outmodel"] = output_models
        like["logfile"]  = log_file
        if force or not os.path.isfile(output_models):
            like.run()
            like.logFileOpen()
        elif os.path.isfile(output_models):
            ml_models = gammalib.GModels(output_models)
        else:
            raise Exception("Cannot proceed with ctlike")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_models)):
            like.save()
            saved = True
        if saved and self.verbosity > 1:
            print("File {} created.".format(output_models), file=sys.stderr)
        return like


