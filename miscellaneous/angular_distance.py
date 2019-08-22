import math as m
import astropy.units as u
from astropy.coordinates import SkyCoord as sc

# coords(degrees)
dat = [ { "ra": 83.6331, "dec": 22.0145 },
        { "ra": 85.6331, "dec": 23.0145 } ]

# https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
star = [ sc(ra=d["ra"], dec=d["dec"], unit="deg", frame="icrs") for d in dat ]
print("s0:", star[0].spherical)
print("s1:", star[1].spherical)

print("Separation [deg]: %.4f" % star[1].separation(star[0]).deg)

for i in range(len(star)):
    print("star", i+1, ":", star[i].to_string("hmsdms"), "||", star[1].ra.hour, star[1].dec.deg, "||", star[1].ra.radian, star[1].dec.radian)

# longitude // theta [0, 2pi]
lon = [ s.ra.radian for s in star ]
# latitude // phi [0, pi]
lat = [ s.dec.radian for s in star ]

# explicit formula: https://en.wikipedia.org/wiki/Great-circle_distance
# phi => declination, lambda => ra
# arccos( sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(delta_lambda) )
# t = m.sin(star[0].dec.radian) * m.sin(star[1].dec.radian) + m.cos(star[0].dec.radian) * m.cos(star[1].dec.radian) * m.cos((star[0].ra.radian-star[1].ra.radian))
t = m.sin(lat[0]) * m.sin(lat[1]) + m.cos(lat[0]) * m.cos(lat[1]) * m.cos(lon[0] - lon[1])
print("Distance [deg]: %.4f" % m.degrees(m.acos(t)))

# Need accuracy about coords (what is theta? what is phi?)
# valentina fioretti code: https://github.com/vfioretti/PythonStuff/blob/master/sph_distance.py
# d_rad = arccos(sin(theta1)sin(theta2)cos(phi1 - phi2) + cos(theta1)cos(theta2)) 
# sph_dist = ((180.0)/math.pi)*math.acos( (math.sin(theta1)*math.sin(theta2)*math.cos(phi1 - phi2) + math.cos(theta1)*math.cos(theta2)) )
phi1, phi2 = lon[0], lon[1]
theta1, theta2 = lat[0], lat[1]
sph_dist = m.sin(theta1) * m.sin(theta2) * m.cos(phi1 - phi2) + m.cos(theta1)*m.cos(theta2)
print("Spherical distance[deg]: %.4f" % m.degrees(m.acos(sph_dist)))
