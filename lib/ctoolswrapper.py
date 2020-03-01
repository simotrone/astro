# Copyright 2019,2020 Simone Tampieri
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import cscripts
import ctools
import gammalib
import os
import sys
import logging
logger = logging.getLogger('ctoolswrapper')

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
            ['nthreads',   10],
        ]
        for f in fields:
            field, default = f
            if default == 'MANDATORY':
                self.__dict__[field] = args[field]
            else:
                self.__dict__[field] = args[field] if field in args and args[field] else default
        logger.setLevel(logging.WARNING - 10*verbosity)

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
        sim["nthreads"] = self.nthreads
        if force or not os.path.isfile(events_file):
            sim.logFileOpen()
            sim.run()
        else:
            container = gammalib.GObservations()
            gcta_obs = gammalib.GCTAObservation(events_file)
            container.append(gcta_obs)
            sim.obs(container)
            sim.obs().models(gammalib.GModels(model_file))
        saved = False
        if (save and force) or (save and not os.path.isfile(events_file)):
            sim.save()
            saved = True
            logger.info("Events file '{}' saved. time [{}-{}]".format(sim["outevents"].value(), time[0], time[1]))
        return sim

    def selection_run(self, input_obs_list, output_obs_list, tmin=0, tmax=None, energy=[None, None], prefix='selected_', log_file='ctselect.log', force=False, save=False):
        if tmin < 0 or tmax < 1 or tmin >= tmax:
            raise Exception("ctselect needs better tmin/tmax")
        select = ctools.ctselect()
        if isinstance(input_obs_list, gammalib.GObservations):
            select.obs(input_obs_list)
        elif os.path.isfile(input_obs_list):
            # observations list from file
            select["inobs"] = input_obs_list
        else:
            raise Exception('Cannot understand input obs list for csphagen')
        select["rad"] = "INDEF"
        select["ra"] = "INDEF"
        select["dec"] = "INDEF"
        select["tmin"] = tmin
        select["tmax"] = tmax
        select["emin"] = energy[0] if energy[0] else "INDEF"
        select["emax"] = energy[1] if energy[1] else "INDEF"
        select["prefix"] = prefix
        select["outobs"] = output_obs_list
        select["logfile"] = log_file
        if force or not os.path.isfile(output_obs_list):
            select.logFileOpen()
            select.run()
        elif os.path.isfile(output_obs_list):
            basename, ext = os.path.splitext(output_obs_list)
            if ext == '.xml':
                select = ctools.ctselect(gammalib.GObservations(output_obs_list))
            else: # .fits
                container = gammalib.GObservations()
                gcta_obs = gammalib.GCTAObservation(output_obs_list)
                container.append(gcta_obs)
                select.obs(container)
        else:
            raise Exception("Cannot proceed with ctselect")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_obs_list)):
            select.save() # why it doesn't save anytime?
            saved = True
            logger.info("Files {} created.".format(output_obs_list))
        return select

    def csphagen_run(self, input_obs_list, input_model=None, source_rad=0.2, energy=[None, None], output_obs_list='onoff_obs.xml', output_model='onoff_result.xml', log_file='csphagen.log', prefix='onoff', ebinalg="LIN", enumbins=1, stacked=False, force=False, save=False):
        phagen = cscripts.csphagen()
        if isinstance(input_obs_list, gammalib.GObservations):
            phagen.obs(input_obs_list)
            if input_model is not None:
                phagen.obs().models(input_model)
        elif os.path.isfile(input_obs_list):
            # observations list from file
            phagen["inobs"] = input_obs_list
            phagen["inmodel"] = input_model
        else:
            raise Exception('Cannot understand input obs list for csphagen')
        phagen["srcname"] = self.name
        phagen["caldb"]   = self.caldb
        phagen["irf"]     = self.irf
        phagen["ebinalg"]  = ebinalg
        phagen["emin"]     = energy[0] if energy[0] else self.energy_min
        phagen["emax"]     = energy[1] if energy[1] else self.energy_max
        phagen["enumbins"] = enumbins
        phagen["coordsys"] = "CEL"
        phagen["ra"]    = self.ra
        phagen["dec"]   = self.dec
        phagen["rad"]   = source_rad
        phagen["stack"] = stacked
        phagen["bkgmethod"] = "REFLECTED"
        phagen["outobs"]   = output_obs_list
        phagen["outmodel"] = output_model
        phagen["prefix"]   = prefix
        phagen["logfile"]  = log_file
        phagen["nthreads"] = self.nthreads
        if force or not os.path.isfile(output_obs_list) or not os.path.isfile(output_model):
            phagen.logFileOpen()
            phagen.run()
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
            logger.info("Files {}, {} created.".format(output_obs_list, output_model))
        return phagen

    def ctlike_run(self, input_obs_list, input_models=None, output_models='ml_result.xml', log_file='ctlike.log', force=False, save=False):
        like = ctools.ctlike()
        if isinstance(input_obs_list, gammalib.GObservations):
            like.obs(input_obs_list)
            if input_models is not None:
                like.obs().models(input_models)
        elif os.path.isfile(input_obs_list) and os.path.isfile(input_models):
            # observations list from file
            like["inobs"] = input_obs_list
            like["inmodel"] = input_models
        else:
            raise Exception('Cannot understand input obs list for ctlike')
        like["outmodel"] = output_models
        like["logfile"]  = log_file
        like["nthreads"] = self.nthreads
        if force or not os.path.isfile(output_models):
            like.logFileOpen()
            like.run()
        elif os.path.isfile(output_models):
            ml_models = gammalib.GModels(output_models)
            like.obs().models(ml_models)
        else:
            raise Exception("Cannot proceed with ctlike")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_models)):
            like.save()
            saved = True
            logger.info("File {} created.".format(output_models))
        return like

    def csspec_run(self, input_obs_list, input_models=None, enumbins=20, output_file='spectrum.fits', log_file='csspec.log', force=False, save=False):
        spec = cscripts.csspec()
        if isinstance(input_obs_list, gammalib.GObservations):
            spec.obs(input_obs_list)
        elif os.path.isfile(input_obs_list) and os.path.isfile(input_models):
            # observations list from file
            spec["inobs"] = input_obs_list
            spec["inmodel"] = input_models
        else:
            raise Exception('Cannot understand input obs list for csspec')
        spec["srcname"] = self.name
        spec["caldb"]   = self.caldb
        spec["irf"]     = self.irf
        spec["method"] = "AUTO"
        spec["emin"] = 0.03
        spec["emax"] = 150.0
        spec['ebinalg']  = "LOG"
        spec["enumbins"] = enumbins
        spec['calc_ts']   = True
        spec['calc_ulim'] = True
        spec['outfile']  = output_file
        spec["logfile"]  = log_file
        spec["nthreads"] = self.nthreads
        if force or not os.path.isfile(output_file):
            spec.logFileOpen()
            spec.run()
        elif os.path.isfile(output_file):
            spec._fits = gammalib.GFits(output_file)
        else:
            raise Exception("Cannot proceed with csspec")
        saved = False
        if (save and force) or (save and not os.path.isfile(output_file)):
            spec.save()
            saved = True
            logger.info("File {} created.".format(output_file))
        return spec

