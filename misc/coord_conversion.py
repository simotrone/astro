from astropy import units as u
from astropy.coordinates import SkyCoord

# coordinates from wikipedia
# frame icrs
c = SkyCoord(ra="05h34m31.94s", dec="22d00m52.2s")
print("ICRS: ", c)

print("   Degree coord: %s" % c.to_string('dms'))
print("Hours+deg coord: %s" % c.to_string('hmsdms'))
print("  Decimal coord: %s" % c.to_string('decimal'))
print("        RA  deg: %.4f" % c.ra.deg)
print("        DEC deg: %.4f" % c.dec.deg)

gal = c.galactic
print("\n")
print(" Galactic: ", gal)
print("   Degree: ", gal.to_string('dms'))
print("Hours+Deg: ", gal.to_string('hmsdms'))
print("  Decimal: ", gal.to_string('decimal'))

fk5 = c.fk5
print("\n")
print("      FK5: ", fk5)
print("   Degree: ", fk5.to_string('dms'))
print("Hours+Deg: ", fk5.to_string('hmsdms'))
print("  Decimal: ", fk5.to_string('decimal'))

