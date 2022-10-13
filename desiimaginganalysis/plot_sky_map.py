from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned
import healpy as hp
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# If False, only plot the footprint
DISPLAY_TARGETS = True

REGION = "north"

plt.rcParams['font.size'] = '12'

nside = 512
npix = hp.nside2npix(nside)

pixel_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)

# Load footprints
with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_north_footprint_for_nside_{nside}.npy", "rb") as f:
    footprint_north = np.load(f)
with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_south_footprint_for_nside_{nside}.npy", "rb") as f:
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

m = np.zeros(npix)

m[footprint_north] = np.nan
m[footprint_south] = np.nan

# Make custom cmap,
# with white as the first color
# so the background values are zero
cmap_jet = get_cmap("jet").copy()
cmap = LinearSegmentedColormap.from_list(
    "",
    ["white"] + list(
        cmap_jet(np.linspace(0, 1, 7))
    ),
)

# If requested, overplot histogram with target counts in each pixel
if DISPLAY_TARGETS:
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_HPXPIXEL_HPXNSIDE_{nside}.npy", "rb") as f:
        target_pixels = np.load(f)

    idx, counts = np.unique(target_pixels, return_counts=True)

    # fill the fullsky map
    m[idx] = counts / pixel_area # Convert from counts to target density (deg^-2)

hp.mollview(
    m,
    fig=0,
    nest=True,
    title=f"Legacy Survey (LS) DR9 Footprint, nside = {nside}",
    unit=r"BGS Target Density [$deg^{-2}$]",
    rot=[120, 0],
    badcolor="gray",
    cbar=True,
    cmap=cmap,
    min=0,
)

# Display lines for Dec=+32 and Dec=-18
# When using projscatter, with lonlat=True,
# first input RA (from 0 to 360), then Dec (from -90 to +90)
ras = np.linspace(0, 360, 700)
north_south_split_radians = 32
desi_southmost_radians = -15

hp.projscatter(ras, [north_south_split_radians]*len(ras), lonlat=True, c="black", s=2, label="N/S Split in LS")
hp.projscatter(ras, [desi_southmost_radians]*len(ras), lonlat=True, c="red", s=2, label="DESI Southmost")


hp.graticule()

plt.legend(
    loc="lower right",
    bbox_to_anchor=(1.02, -0.08)
)

plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/footprint.pdf")
plt.show()
