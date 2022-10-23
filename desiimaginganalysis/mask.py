import healpy as hp
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from desitarget.io import desitarget_resolve_dec
from desitarget.geomask import nside2nside

def mask(pixel_ids, nside, region="both"):
    # Given an array of pixel IDs at some NSIDE (in NESTED ordering),
    # return an array of True and False values indicating whether
    # each pixel is within the mask.
    # Can choose north, south, or both
    # for the mask region used.
    if region not in {"north", "south", "both"}:
        raise ValueError("region must be north, south, or both.")

    with open ("/cluster/scratch/lmachado/DataProducts/masks/north_mask.npy", "rb") as f:
        north_mask = np.load(f)
    with open ("/cluster/scratch/lmachado/DataProducts/masks/south_mask.npy", "rb") as f:
        south_mask = np.load(f)

    NSIDE_MASKS = hp.get_nside(north_mask)

    # Change mask maps to the nside of the pixel_ids passed into the function
    north_mask = hp.pixelfunc.ud_grade(
        north_mask,
        nside_out=nside,
        order_in="NEST",
        order_out="NEST",
        dtype=int,
    )
    south_mask = hp.pixelfunc.ud_grade(
        south_mask,
        nside_out=nside,
        order_in="NEST",
        order_out="NEST",
        dtype=int,
    )

    north_results = north_mask[pixel_ids]
    south_results = south_mask[pixel_ids]

    if region == "north":
        return north_results
    if region == "south":
        return south_results
    if region == "both":
        return north_results | south_results

# Running this file directly will
# execute the mask creation code, which
# plots a few validation figures, and
# saves the resulting masks into npy files.
if __name__ == "__main__":
    # Load footprint and targets
    NSIDE_FILES = 512
    NSIDE = 64

    NPIX_FILES = hp.nside2npix(NSIDE_FILES)
    NPIX = hp.nside2npix(NSIDE)

    NORTH_FOOTPRINT_FILE = f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_north_footprint_for_nside_{NSIDE_FILES}.npy"
    SOUTH_FOOTPRINT_FILE = f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_south_footprint_for_nside_{NSIDE_FILES}.npy"

    NORTH_TARGET_PIXELS_FILE = f"/cluster/scratch/lmachado/DataProducts/targets/north/targets_HPXPIXEL_HPXNSIDE_{NSIDE_FILES}.npy"
    SOUTH_TARGET_PIXELS_FILE = f"/cluster/scratch/lmachado/DataProducts/targets/south/targets_HPXPIXEL_HPXNSIDE_{NSIDE_FILES}.npy"

    with open(NORTH_FOOTPRINT_FILE, "rb") as f:
        north_footprint_pixels = np.load(f)
    with open(NORTH_TARGET_PIXELS_FILE, "rb") as f:
        north_target_pixels = np.load(f)
    with open(SOUTH_FOOTPRINT_FILE, "rb") as f:
        south_footprint_pixels = np.load(f)
    with open(SOUTH_TARGET_PIXELS_FILE, "rb") as f:
        south_target_pixels = np.load(f)

    full_footprint_pixels = np.concatenate((north_footprint_pixels, south_footprint_pixels))
    full_target_pixels = np.concatenate((north_target_pixels, south_target_pixels))

    north_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, north_footprint_pixels)
    north_target_pixels = nside2nside(NSIDE_FILES, NSIDE, north_target_pixels)
    south_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, south_footprint_pixels)
    south_target_pixels = nside2nside(NSIDE_FILES, NSIDE, south_target_pixels)

    full_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, full_footprint_pixels)
    full_target_pixels = nside2nside(NSIDE_FILES, NSIDE, full_target_pixels)

    # Obtain pixels in correct region
    # Convert cutoff from degrees to colatitude
    north_south_cutoff = np.pi/2 - np.radians(desitarget_resolve_dec())

    # North
    ring_north_pixels = hp.query_strip(
        NSIDE,
        theta1=0,
        theta2=north_south_cutoff,
        inclusive=True,
        nest=False,
    )

    ring_north_map = np.zeros(NPIX)
    ring_north_map[ring_north_pixels] = 1

    nest_north_map = hp.reorder(ring_north_map, r2n=True)

    nest_north_pixels = np.where(nest_north_map == 1)[0]

    # South
    ring_south_pixels = hp.query_strip(
        NSIDE,
        theta1=north_south_cutoff,
        theta2=np.pi,
        inclusive=True,
        nest=False,
    )

    ring_south_map = np.zeros(NPIX)
    ring_south_map[ring_south_pixels] = 1

    nest_south_map = hp.reorder(ring_south_map, r2n=True)

    nest_south_pixels = np.where(nest_south_map == 1)[0]

    # Sanity check: do these region maps add up to the whole sky?
    assert NPIX == len(set(ring_north_pixels).union(set(ring_south_pixels)))
    assert NPIX == len(set(nest_north_pixels).union(set(nest_south_pixels)))

    # Remove pixels from footprint that do not belong to the desired region
    north_footprint_pixels = np.array(list(
        set(north_footprint_pixels).intersection(set(nest_north_pixels))
    ))
    south_footprint_pixels = np.array(list(
        set(south_footprint_pixels).intersection(set(nest_south_pixels))
    ))

    # Sanity check: do north and south footprints, when combined,
    # touch nicely? Or is there an undesired gap between them?
    m = np.zeros(NPIX)
    m[north_footprint_pixels] = 1
    m[south_footprint_pixels] = 2
    hp.mollview(m, nest=True, rot=[120, 0], title="Footprints")
    hp.graticule()
    plt.show()


    ###############
    # Combined mask
    ###############
    # Get all neighbor pixels of all target pixels
    # Note: since phi=None, theta is taken as a list of pixel IDs
    all_neighbor_pixels = hp.pixelfunc.get_all_neighbours(
        NSIDE,
        theta=full_target_pixels,
        phi=None,
        nest=True,
    )

    all_neighbor_pixels = all_neighbor_pixels.flatten()

    # Remove target pixels from the neighbor list
    right_outside_pixels = np.array(list(
        set(all_neighbor_pixels).difference(set(full_target_pixels))
    ))

    # Obtain neighbors of the pixels right outside
    right_outside_neighbors = hp.pixelfunc.get_all_neighbours(
        NSIDE,
        theta=right_outside_pixels,
        phi=None,
        nest=True,
    )

    right_outside_neighbors = right_outside_neighbors.flatten()

    # Remove such neighbors from the mask,
    # since they are in the boundary
    pixels_in_mask = set(full_target_pixels).difference(set(right_outside_neighbors))

    pixels_in_mask = np.array(list(pixels_in_mask))

    # Visualize mask and targets together at NSIDE
    # Make custom cmap,
    # with white as the first color
    # so the background values are zero
    cmap = LinearSegmentedColormap.from_list(
        "",
        ["white", "blue", "gray"]
    )

    m = np.zeros(NPIX)
    m[full_target_pixels] = 1
    m[pixels_in_mask] = 2

    hp.mollview(m, nest=True, rot=[120, 0], title="Mask and target area", cbar=False, cmap=cmap)

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
    plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/mask.pdf")
    plt.show()

    # Sanity check:
    # Is the mask totally inside the footprint?
    mask_footprint_intersection = set(pixels_in_mask).intersection(set(full_footprint_pixels))
    # If the lenghts are not the same, the mask is not entirely contained in the footprint
    assert len(mask_footprint_intersection) == len(set(pixels_in_mask))

    # Sanity check:
    # Is the mask totally inside the targets area?
    mask_target_intersection = set(pixels_in_mask).intersection(set(full_target_pixels))
    # If the lenghts are not the same, the mask is not entirely contained in the target area
    assert len(mask_target_intersection) == len(set(pixels_in_mask))

    # Sanity check:
    # Are there any holes in the gray footprint?
    footprint_holes = set(pixels_in_mask).difference(set(full_footprint_pixels))
    assert len(footprint_holes) == 0

    # Sanity check:
    # Are there any holes in the targets pixels?
    target_holes = set(pixels_in_mask).difference(set(full_target_pixels))
    assert len(target_holes) == 0

    m = np.zeros(NPIX)
    m[pixels_in_mask] = 1
    hp.mollview(m, nest=True, rot=[120, 0], title="Mask", cbar=False, cmap=cmap)
    hp.graticule()
    plt.show()

    # Compute area percentages
    footprint_loss = (len(set(full_footprint_pixels)) - len(set(pixels_in_mask))) / len(set(full_footprint_pixels))
    target_loss = (len(set(full_target_pixels)) - len(set(pixels_in_mask))) / len(set(full_target_pixels))

    print("Gray footprint loss:", footprint_loss)
    print("Target area loss:", target_loss)

    # Obtain northern and southern masks by taking intersection with each region
    pixels_in_north_mask = np.array(list(
        set(pixels_in_mask).intersection(set(nest_north_pixels))
    ))
    pixels_in_south_mask = np.array(list(
        set(pixels_in_mask).intersection(set(nest_south_pixels))
    ))

    # Sanity check: when plotting both masks together,
    # do they touch nicely at the border?
    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_north_mask] = 1
    m[pixels_in_south_mask] = 2
    hp.mollview(m, nest=True, rot=[120, 0], title="North and South masks together", cmap=cmap, cbar=False)
    hp.graticule()
    plt.show()

    # Finally, save the resulting masks as bit arrays
    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_north_mask] = 1
    with open ("/cluster/scratch/lmachado/DataProducts/masks/north_mask.npy", "wb") as f:
        np.save(f, m)

    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_south_mask] = 1
    with open ("/cluster/scratch/lmachado/DataProducts/masks/south_mask.npy", "wb") as f:
        np.save(f, m)
