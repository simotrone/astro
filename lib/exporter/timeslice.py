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

from astropy.io import fits
from defusedxml.ElementTree import parse
from lib.exporter.csv import CSVExporter as csvex
import os
import sys

class TimeSliceExporter:
    # Reference: http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#file-function
    def __init__(self, input_filename, template_model=None, savings_dir=None, verbosity=0, tmax=None):
        self.input_filename = input_filename
        self.hdul = fits.open(input_filename)
        self.verbosity = verbosity
        if self.verbosity > 0:
            self.hdul.info()
        if tmax is not None and tmax < 1:
            raise Exception('Need a non-zero positive tmax')
        self.tmax = tmax

        self.model_filename = template_model
        self.model_tree = None
        if self.model_filename:
            self.model_tree = parse(self.model_filename)

        self.savings_dir = '.'
        if savings_dir:
            if os.path.isabs(savings_dir):
                self.savings_dir = savings_dir
            else:
                self.savings_dir = os.path.join(self.savings_dir, savings_dir)
            try:
                os.makedirs(self.savings_dir)
                if self.verbosity > 0:
                    print("Created the data dir: {}".format(self.savings_dir), file=sys.stderr)
            except FileExistsError as e:
                if self.verbosity > 1:
                    print("The data dir already exists", file=sys.stderr)

    # times in secs
    # energies in GeV
    # spectra in ph/cm2/s/GeV
    # Important: we need to convert GeV in MeV (see ref.)
    def export(self, force=False, ebl=False):
        times    = self.hdul['TIMES'].data
        energies = self.hdul['ENERGIES'].data
        spectra  = None
        if ebl:
            spectra = self.hdul['EBL GILMORE'].data
        else:
            spectra = self.hdul['SPECTRA'].data
        mev_energies = [ene[0]*1000 for ene in energies]
        done = []
        for i, tsec in enumerate(times):
            # writing energies/spectra tsv
            mev_spectra = [f/1000 for f in spectra[i]]
            time_slice_filename = os.path.join(self.savings_dir, "spec_tbin{0:02d}.tsv".format(i))
            if force or not os.path.isfile(time_slice_filename):
                self.write_energies_spectra_tsv(time_slice_filename, mev_energies, mev_spectra)
            elif self.verbosity > 1:
                print("The slice file {} already exists".format(time_slice_filename), file=sys.stderr)
            # writing model xml if template was provided
            xml_slice_filename = None
            if self.model_filename:
                base, ext = os.path.splitext(os.path.basename(self.model_filename))
                xml_slice_filename = os.path.join(self.savings_dir, base+'_tbin{0:02d}'.format(i)+ext)
                if force or not os.path.isfile(xml_slice_filename):
                    self.write_linked_model(self.model_tree, os.path.basename(time_slice_filename), xml_slice_filename)
                elif self.verbosity > 1:
                    print("The xml file {} already exists".format(xml_slice_filename), file=sys.stderr)
            done.append({
                "id": i,
                "tsec": tsec[0],
                "ene_flux_file": time_slice_filename,
                "model_file": xml_slice_filename,
            })
            if self.verbosity > 1:
                print('slice {0:2d} {1:15f} sec > {2}{3}'.format(i, tsec[0], time_slice_filename, ', '+xml_slice_filename if xml_slice_filename else ''), file=sys.stderr)
            # test this at the end because we want  time slot more over the tmax
            if self.tmax is not None and self.tmax < tsec[0]:
                break
        return done

    # http://cta.irap.omp.eu/gammalib-devel/doxygen/classGModelSpectralFunc.html#a922f555636e58e519871913bcd76c5d8
    # http://cta.irap.omp.eu/gammalib-devel/doxygen/classGCsv.html#details
    # GCsv: This class implements a table of std::string elements that is loaded
    #       from a comma-separated value ASCII file. The comma-separation string
    #       can be specified upon loading of the file (by default the class
    #       assumes that elements are separated by a white space).
    @classmethod
    def write_energies_spectra_tsv(cls, output_filename, energies, spectra):
        if not len(energies) == len(spectra):
            raise Exception('Need the same number of elements between energies and spectra')
        data = [ [energy, spectra[i]] for i, energy in enumerate(energies) ]
        csvex.save(output_filename, data, delimiter=" ")

    @classmethod
    def write_linked_model(cls, model_tree, ref_file, output_xml_fn):
        if not model_tree:
            raise Exception('Need a defined xml model tree')
        source_library = model_tree.getroot()
        for source in source_library:
            spectrum = source.find('spectrum')
            if spectrum.attrib.get('type') != 'FileFunction':
                continue
            spectrum.set('file', ref_file)
            break # we can stop at first occurrency
        model_tree.write(output_xml_fn, xml_declaration=True, encoding="UTF-8")

