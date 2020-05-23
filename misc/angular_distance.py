import math as m
import astropy.units as u
from astropy.coordinates import SkyCoord as sc

data = [
    # source                           ,  pointing
    # Crab
    #[ { "ra": 83.6331, "dec": 22.0145 }, { "ra": 83.6331, "dec": 22.5145 } ], # +0.5 dec
    #[ { "ra": 83.6331, "dec": 22.0145 }, { "ra": 84.1331, "dec": 22.0145 } ], # +0.5 ra
    #[ { "ra": 83.6331, "dec": 22.0145 }, { "ra": 85.6331, "dec": 23.0145 } ], # +2.0 ra, +1 dec
    # GRB
    [ { "ra": 33.057, "dec": -51.841 }, { "ra": 33.857,  "dec": -51.841 } ], # +0.8 ra
    [ { "ra": 33.057, "dec": -51.841 }, { "ra": 32.257, 'dec': -51.841  } ],
    [ { 'ra': 32.253, 'dec': -52.335 }, { "ra": 32.257, 'dec': -51.841  } ],
    [ { 'ra': 31.457, 'dec': -51.836 }, { "ra": 32.257, 'dec': -51.841  } ],
    [ { 'ra': 32.261, 'dec': -51.347 }, { "ra": 32.257, 'dec': -51.841  } ],

    # p0_R, +0.5 ra, distanza 0.3deg
    #[ { "ra": 33.057,  "dec": -51.841 }, { "ra": 33.057,  "dec": -51.341 } ], # +0.5 dec
    # [ { "ra": 33.057,  "dec": -51.841 }, { "ra": 33.557,  "dec": -51.841 } ], # +0.5 ra
    # [ { "ra": 33.555,  "dec": -51.532 }, { "ra": 33.557,  "dec": -51.841 } ],
    # [ { "ra": 34.056,  "dec": -51.838 }, { "ra": 33.557,  "dec": -51.841 } ],
    # [ { "ra": 33.558,  "dec": -52.149 }, { "ra": 33.557,  "dec": -51.841 } ],

    # p1_U, +0.5 dec, distanza 0.5deg
    #[ { "ra": 33.057,  "dec": -51.841 }, { "ra": 33.057,  "dec": -52.341 } ], # +0.5 dec
    #[ { "ra": 33.875,  "dec": -52.338 }, { "ra": 33.057,  "dec": -52.341 } ],
    #[ { "ra": 33.056,  "dec": -52.841 }, { "ra": 33.057,  "dec": -52.341 } ],
    #[ { "ra": 32.238,  "dec": -52.338 }, { "ra": 33.057,  "dec": -52.341 } ],


    # ctools csphagen simulation
    # [ { "ra": 33.057,  "dec": -51.841 }, { "ra": 33.857,  "dec": -51.841 } ], # +0.8 ra # source, pnt
    # [ { "ra": 34.029,  "dec": -51.359 }, { "ra": 33.857,  "dec": -51.841 } ],
    # [ { "ra": 34.572,  "dec": -51.622 }, { "ra": 33.857,  "dec": -51.841 } ],
    # [ { "ra": 34.583,  "dec": -52.051 }, { "ra": 33.857,  "dec": -51.841 } ],
    # [ { "ra": 34.041,  "dec": -52.322 }, { "ra": 33.857,  "dec": -51.841 } ],

    # [ { "ra": 33.057,  "dec": -51.841 }, { "ra": 34.057,  "dec": -51.841 } ], # +1.0 ra # source, pnt
    # [ { "ra": 33.878,  "dec": -51.233 }, { "ra": 34.057,  "dec": -51.841 } ],
    # [ { "ra": 34.545,  "dec": -51.302 }, { "ra": 34.057,  "dec": -51.841 } ],
    # [ { "ra": 34.989,  "dec": -51.622 }, { "ra": 34.057,  "dec": -51.841 } ],
    # [ { "ra": 35.003,  "dec": -52.044 }, { "ra": 34.057,  "dec": -51.841 } ],
    # [ { "ra": 34.568,  "dec": -52.372 }, { "ra": 34.057,  "dec": -51.841 } ],
    # [ { "ra": 33.887,  "dec": -52.450 }, { "ra": 34.057,  "dec": -51.841 } ],

    # distance variability with RA fixed, declination move
    # [ { "ra": 0.0,  "dec": 0.0 }, { "ra": 0.5,  "dec": 0.0 } ],
    # [ { "ra": 0.0,  "dec": 10.0 }, { "ra": 0.5,  "dec": 10.0 } ],
    # [ { "ra": 0.0,  "dec": 20.0 }, { "ra": 0.5,  "dec": 20.0 } ],
    # [ { "ra": 0.0,  "dec": 30.0 }, { "ra": 0.5,  "dec": 30.0 } ],
    # [ { "ra": 0.0,  "dec": 40.0 }, { "ra": 0.5,  "dec": 40.0 } ],
    # [ { "ra": 0.0,  "dec": 50.0 }, { "ra": 0.5,  "dec": 50.0 } ],
    # [ { "ra": 0.0,  "dec": -50.0 }, { "ra": 0.5,  "dec": -50.0 } ],
    # [ { "ra": 0.0,  "dec": 60.0 }, { "ra": 0.5,  "dec": 60.0 } ],
    # [ { "ra": 0.0,  "dec": 70.0 }, { "ra": 0.5,  "dec": 70.0 } ],
    # [ { "ra": 0.0,  "dec": 80.0 }, { "ra": 0.5,  "dec": 80.0 } ],
]

for dat in data:
    # https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    star = [ sc(ra=d["ra"], dec=d["dec"], unit="deg", frame="icrs") for d in dat ]
    print("s0:", star[0].spherical)
    print("s1:", star[1].spherical)
    print("  Separation [deg]: %.4f" % star[1].separation(star[0]).deg)

    for i in range(len(star)):
        print("  star", i+1, ":", star[i].to_string("hmsdms"), "||", star[1].ra.hour, star[1].dec.deg, "||", star[1].ra.radian, star[1].dec.radian)

    # longitude // theta [0, 2pi]
    lon = [ s.ra.radian for s in star ]
    # latitude // phi [0, pi]
    lat = [ s.dec.radian for s in star ]

    # explicit formula: https://en.wikipedia.org/wiki/Great-circle_distance
    # phi => declination, lambda => ra
    # arccos( sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(delta_lambda) )
    # t = m.sin(star[0].dec.radian) * m.sin(star[1].dec.radian) + m.cos(star[0].dec.radian) * m.cos(star[1].dec.radian) * m.cos((star[0].ra.radian-star[1].ra.radian))
    t = m.sin(lat[0]) * m.sin(lat[1]) + m.cos(lat[0]) * m.cos(lat[1]) * m.cos(lon[0] - lon[1])
    ds = m.acos(t)
    print("  Distance [deg]: %.4f" % m.degrees(ds))

    if False:
        # Need accuracy about coords (what is theta? what is phi?)
        # valentina fioretti code: https://github.com/vfioretti/PythonStuff/blob/master/sph_distance.py
        # d_rad = arccos(sin(theta1)sin(theta2)cos(phi1 - phi2) + cos(theta1)cos(theta2)) 
        # sph_dist = ((180.0)/math.pi)*math.acos( (math.sin(theta1)*math.sin(theta2)*math.cos(phi1 - phi2) + math.cos(theta1)*math.cos(theta2)) )
        phi1, phi2 = lon[0], lon[1]
        theta1, theta2 = lat[0], lat[1]
        sph_dist = m.sin(theta1) * m.sin(theta2) * m.cos(phi1 - phi2) + m.cos(theta1)*m.cos(theta2)
        sph_dist_arcos = m.acos(sph_dist)
        print("  Spherical distance[deg]: %.4f" % m.degrees(sph_dist_arcos))
