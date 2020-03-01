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

import csv
import os

class CSVExporter:
    # {'id': 0, 'tsec': 0.1, 'ene_flux_file': './data/spec_00.tsv', 'model_file': './data/run0406_ID000126_00.xml'},
    # {'id': 1, 'tsec': 0.12589253, 'ene_flux_file': './data/spec_01.tsv', 'model_file': './data/run0406_ID000126_01.xml'},
    # {'id': 2, 'tsec': 0.15848932, 'ene_flux_file': './data/spec_02.tsv', 'model_file': './data/run0406_ID000126_02.xml'}, 
    # ...
    @staticmethod
    def save(output_filename, data, headers=None, delimiter=" "):
        if not isinstance(data, list):
            raise Exception('Need a data list to write')
        with open(output_filename, mode='w', newline="\n") as fh:
            writer = None
            csv_options = { 'delimiter': delimiter, 'quotechar':'"', 'quoting':csv.QUOTE_MINIMAL }
            if headers:
                writer = csv.DictWriter(fh, fieldnames=headers, **csv_options)
                writer.writeheader()
            else:
                writer = csv.writer(fh, **csv_options)
            writer.writerows(data)

    # @staticmethod
    def load(input_filename, header=False, sep=' '):
        if not os.path.isfile(input_filename):
            raise Exception('Cannot read {}'.format(input_filename))
        data = []
        with open(input_filename, newline='\n') as fh:
            reader = csv.reader(fh, delimiter=sep)

            if header:
                headers_row = next(reader, None)

            for row in reader:
                if header:
                    d = {}
                    for i, h in enumerate(headers_row):
                        d[h] = row[i]
                else:
                    d = row
                data.append(d)
        return data
