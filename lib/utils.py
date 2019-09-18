import csv
import math

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

def read_timeslices_tsv(filename):
    ts = []
    with open(filename, mode='r', newline='\n') as fh:
        reader = csv.reader(fh, delimiter='\t')
        headers = next(reader)
        for row in reader:
            ts.append(dict(zip(headers, row)))
    return ts
