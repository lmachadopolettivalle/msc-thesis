from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned
import healpy as hp
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np

plt.rcParams['font.size'] = '12'

nside = 512
npix = hp.nside2npix(nside)

# Load footprints
with open(f"pixels_in_north_footprint_for_nside_{nside}.npy", "rb") as f:
    footprint_north = np.load(f)
with open(f"pixels_in_south_footprint_for_nside_{nside}.npy", "rb") as f:
    footprint_south = np.load(f)

# Compute areas in square degrees of each mask
pixel_area = hp.nside2pixarea(nside, degrees=True)
area_footprint_north = pixel_area * len(footprint_north)
area_footprint_south = pixel_area * len(footprint_south)

footprint_intersection = np.array(list(set(footprint_north).intersection(set(footprint_south))))
area_footprint_intersection = pixel_area * len(footprint_intersection)

print("Areas in square degrees:")
print("north: ", area_footprint_north)
print("south: ", area_footprint_south)
print("intersection: ", area_footprint_intersection)

# Make empty mollweide projection plot, with equatorial coordinates
#ax = init_sky(galactic_plane_color=None, ecliptic_plane_color=None)


m = np.zeros(npix)

m[footprint_north] = 0.5
m[footprint_south] = 0.5

cmap = get_cmap("gray_r").copy()

hp.mollview(
    m,
    nest=True,
    title=f"Legacy Survey DR9 Footprint, nside = {nside}",
    rot=[120,0],
    cbar=False,
    cmap=cmap,
    min=0,
    max=1,
)

# Display lines for Dec=+32 and Dec=-18
# When using projscatter, with lonlat=True,
# first input RA (from -180 to +180), then Dec (from -90 to +90)
ras = np.linspace(-180, 180, 700)
ras = np.linspace(0, 360, 700)
north_south_split_radians = 32
desi_southmost_radians = -15

hp.projscatter(ras, [north_south_split_radians]*len(ras), lonlat=True, c="black", s=2, label="North-South Split in Legacy Survey")
hp.projscatter(ras, [desi_southmost_radians]*len(ras), lonlat=True, c="red", s=2, label="DESI Southmost Coverage")


hp.graticule()

plt.legend(loc="lower right", bbox_to_anchor=(1, -0.16))

plt.savefig("footprint.pdf")
plt.show()
