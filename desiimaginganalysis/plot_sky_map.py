from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned
import healpy as hp
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

DISPLAY_FOOTPRINTS = {
    "north": False,
    "south": False,
}
DISPLAY_TARGETS = {
    "north": True,
    "south": True,
}

plt.rcParams['font.size'] = '12'

nside = 512
npix = hp.nside2npix(nside)

pixel_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)


def print_footprint_areas():
    # Compute areas in square degrees of each mask
    with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_north_footprint_for_nside_{nside}.npy", "rb") as f:
        footprint_north = np.load(f)
    with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_south_footprint_for_nside_{nside}.npy", "rb") as f:
        footprint_south = np.load(f)

    area_footprint_north = pixel_area * len(footprint_north)
    area_footprint_south = pixel_area * len(footprint_south)

    footprint_intersection = np.array(list(set(footprint_north).intersection(set(footprint_south))))
    area_footprint_intersection = pixel_area * len(footprint_intersection)

    print("Areas in square degrees:")
    print("north: ", area_footprint_north)
    print("south: ", area_footprint_south)
    print("intersection: ", area_footprint_intersection)


# Load footprints
m = np.zeros(npix)

for region, should_display in DISPLAY_FOOTPRINTS.items():
    if should_display:
        with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_{region}_footprint_for_nside_{nside}.npy", "rb") as f:
            footprint = np.load(f)

        m[footprint] = np.nan


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

unit_label = None
# If requested, overplot histogram with target counts in each pixel
for region, should_display in DISPLAY_TARGETS.items():
    if should_display:
        unit_label = r"BGS Target Density [$deg^{-2}$]"

        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_HPXPIXEL_HPXNSIDE_{nside}.npy", "rb") as f:
            target_pixels = np.load(f)

        idx, counts = np.unique(target_pixels, return_counts=True)

        # fill the fullsky map
        m[idx] = counts / pixel_area # Convert from counts to target density (deg^-2)

# Determine if should show colorbar,
# i.e. if will show any targets
DISPLAY_COLORBAR = False
for should_display in DISPLAY_TARGETS.values():
    if should_display:
        DISPLAY_COLORBAR = True
        break

hp.mollview(
    m,
    fig=0,
    nest=True,
    title=f"Legacy Survey (LS) DR9 Footprint, nside = {nside}",
    unit=unit_label,
    rot=[120, 0],
    badcolor="gray",
    cbar=DISPLAY_COLORBAR,
    cmap=cmap,
    min=0,
    max=(None if DISPLAY_COLORBAR else 1), # If no targets, make max=1 to make background white
)

# Display lines for Dec=+32 and Dec=-18
# When using projscatter, with lonlat=True,
# first input RA (from 0 to 360), then Dec (from -90 to +90)
ras = np.linspace(0, 360, 700)
north_south_split_degrees = 32
desi_southmost_degrees = -15

hp.projscatter(ras, [north_south_split_degrees]*len(ras), lonlat=True, c="black", s=2, label="N/S Split in LS")
hp.projscatter(ras, [desi_southmost_degrees]*len(ras), lonlat=True, c="red", s=2, label="DESI Southmost")


hp.graticule()

plt.legend(
    loc="lower right",
    bbox_to_anchor=(1.02, -0.08)
)

plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/targets_and_footprint.pdf")
plt.show()
