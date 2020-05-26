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
import math
from astropy.coordinates import SkyCoord, Angle

def li_ma (n_on, n_off, alpha):
    if n_on <= 0 or n_off <= 0 or alpha == 0:
        return None
    fc = 1 / alpha * (1 + alpha)
    fb = n_on / (n_on + n_off)
    f  = fc * fb
    gc = 1 + alpha
    gb = n_off / (n_on + n_off)
    g  = gc * gb
    first  = n_on * math.log(f)
    second = n_off * math.log(g)
    fullb   = first + second
    return math.sqrt(2) * math.sqrt(fullb)

def read_timeslices_tsv(filename):
    ts = []
    with open(filename, mode='r', newline='\n') as fh:
        reader = csv.reader(fh, delimiter='\t')
        headers = next(reader)
        for row in reader:
            ts.append(dict(zip(headers, row)))
    return ts

def get_angle(input_angle):
    ang = None
    if isinstance(input_angle, Angle):
        ang = input_angle
    elif isinstance(input_angle, float):
        ang = Angle(input_angle, unit='deg')
    else:
        raise Exception('The input parameter must be an Angle or a float for decimal degree.')
    return ang

def get_skycoord(pnt_coord):
    coord = None
    if isinstance(pnt_coord, SkyCoord):
        coord = pnt_coord
    elif isinstance(pnt_coord, dict) and 'ra' in pnt_coord and 'dec' in pnt_coord:
        coord = SkyCoord(ra=pnt_coord['ra'], dec=pnt_coord['dec'], unit='deg', frame='icrs')
    else:
        raise Exception('The input parameter must be a SkyCoord or a { "ra": 12.3, "dec": 45.6 } dictionary.')
    return coord

