# testing coords and projection
import sys
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io  import fits

# input: a skymap
hdu = fits.open(sys.argv[1])[0]
wcs = WCS(hdu.header)
print(wcs)
# distortion
if False:
    wcs.wcs.crval = [73, 89]
    # wcs.wcs.crpix = [0, 0]
    print(wcs)
# print(wcs.to_header())

fig = plt.figure()
ax = fig.add_subplot(projection=wcs)
ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
ax.grid(color='white', ls='solid')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.coords[0].set_major_formatter('d.ddd')
# overlay = ax.get_coords_overlay('fk5')
# overlay = ax.get_coords_overlay('icrs')
# overlay.grid(color='green', ls='dotted')
# overlay[0].set_axislabel('RA (J2000)')
# overlay[1].set_axislabel('Dec (J2000)')

# galactic overlay
overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='cyan', ls='solid')
overlay[0].set_axislabel('Gal Lon')
overlay[1].set_axislabel('Gal Lat')

ax.scatter(200.5, 200.5, marker='+', color='r')
ax.scatter(83.6, 22.0, marker='*', color='k', transform=ax.get_transform('world'))

plt.show()
exit(1)

