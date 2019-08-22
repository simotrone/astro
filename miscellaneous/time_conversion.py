import itertools
from astropy.time import Time

times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00',
         '2019-01-01T12:00:00', '2019-01-01T12:30:00']

t = Time(times, format='isot', scale='utc')
m = Time([51544.5, 51544.511574], format='mjd')
print("t: ", t)
print("m: ", m)

for i in itertools.chain(t, m):
    print("time: %23s, ISO: %23s, JD: %.5f, MJD: %.5f" % (i, i.iso, i.jd, i.mjd))
