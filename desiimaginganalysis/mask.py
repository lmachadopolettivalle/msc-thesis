import healpy as hp
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from desitarget.io import desitarget_resolve_dec
from desitarget.geomask import hp_in_box, nside2nside

from constants import BASS_MzLS, DECaLS_NGC, DECaLS_SGC, ALL_REGIONS

def mask(pixel_ids, nside, regions=ALL_REGIONS):
    # Given an array of pixel IDs at some NSIDE (in NESTED ordering),
    # return an array of True and False values indicating whether
    # each pixel is within the mask.
    # Can choose north, south, or both
    # for the mask region used.
    for region in regions:
        if region not in ALL_REGIONS:
            raise ValueError("Region is not valid.")

    with open ("/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy", "rb") as f:
        BASS_MzLS_mask = np.load(f)
    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_NGC_mask.npy", "rb") as f:
        DECaLS_NGC_mask = np.load(f)
    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_SGC_mask.npy", "rb") as f:
        DECaLS_SGC_mask = np.load(f)

    NSIDE_MASKS = hp.get_nside(BASS_MzLS_mask)

    # Change mask maps to the nside of the pixel_ids passed into the function
    BASS_MzLS_mask = hp.pixelfunc.ud_grade(
        BASS_MzLS_mask,
        nside_out=nside,
        order_in="NEST",
        order_out="NEST",
        dtype=int,
    )
    DECaLS_NGC_mask = hp.pixelfunc.ud_grade(
        DECaLS_NGC_mask,
        nside_out=nside,
        order_in="NEST",
        order_out="NEST",
        dtype=int,
    )
    DECaLS_SGC_mask = hp.pixelfunc.ud_grade(
        DECaLS_SGC_mask,
        nside_out=nside,
        order_in="NEST",
        order_out="NEST",
        dtype=int,
    )

    BASS_MzLS_results = BASS_MzLS_mask[pixel_ids]
    DECaLS_NGC_results = DECaLS_NGC_mask[pixel_ids]
    DECaLS_SGC_results = DECaLS_SGC_mask[pixel_ids]

    results = np.zeros(len(pixel_ids), dtype=int)

    if BASS_MzLS in regions:
        results = results | BASS_MzLS_results
    if DECaLS_NGC in regions:
        results = results | DECaLS_NGC_results
    if DECaLS_SGC in regions:
        results = results | DECaLS_SGC_results

    return results

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

    # Print number of pixels in each area, AFTER downgrading pixels
    print("Number of pixels in north target set, before masking:", len(set(north_target_pixels)))
    print("Number of pixels in south target set, before masking:", len(set(south_target_pixels)))

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

    # Remove small, disjoint patches from mask to make it simpler
    # Regions to be removed: chosen by eye based on initial version of mask
    pixels_to_be_removed = set()
    for radecbox in [
        [210, 240, -30, -10],
        [120, 180, -40, -10.2],
    ]:
        pixels_to_be_removed.update(list(hp_in_box(
            NSIDE,
            radecbox,
            inclusive=False,
        )))

    pixels_in_mask = pixels_in_mask.difference(pixels_to_be_removed)

    # Convert to np.array for visualization
    pixels_in_mask = np.array(list(pixels_in_mask))


    # Visualize mask and targets together at NSIDE
    # Make custom cmap,
    # with white as the first color
    # so the background values are zero
    cmap = LinearSegmentedColormap.from_list(
        "",
        ["white", "blue", "orange", "green", "gray"]
    )

    m = np.zeros(NPIX)
    m[full_target_pixels] = 1
    m[pixels_in_mask] = 4

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

    # Obtain masks for each of the three regions from the north and south masks

    # BASS/MzLS
    # First, intersect full mask with north pixels.
    # Then, intersect result with RA/DEC box.
    BASS_MzLS_box = set(list(hp_in_box(
        NSIDE,
        [60, 210, desitarget_resolve_dec()-1, 90-1e-3],
        inclusive=True,
    ))).union(set(list(
        hp_in_box(
            NSIDE,
            [190, 310, desitarget_resolve_dec()-1, 90-1e-3],
            inclusive=True,
    ))))
    pixels_in_BASS_MzLS = set(pixels_in_mask).intersection(set(nest_north_pixels)).intersection(set(BASS_MzLS_box))

    # DECaLS NGC
    # First, intersect with south pixels.
    # Then, intersect result with RA/DEC box.
    DECaLS_NGC_box = set(list(hp_in_box(
        NSIDE,
        [90, 210, -15, desitarget_resolve_dec()-1],
        inclusive=True,
    ))).union(set(list(
        hp_in_box(
            NSIDE,
            [190, 290, -15, desitarget_resolve_dec()-1],
            inclusive=True,
    ))))

    pixels_in_DECaLS_NGC = set(pixels_in_mask).intersection(set(nest_south_pixels)).intersection(set(DECaLS_NGC_box))

    # DECaLS SGC
    # This is the full mask, minus the BASS/MzLS pixels,
    # minus the DECaLS NGC pixels.
    pixels_in_DECaLS_SGC = set(pixels_in_mask).difference(pixels_in_BASS_MzLS).difference(pixels_in_DECaLS_NGC)

    # Sanity check: do the 3 masks have all the pixels
    # that used to be in the original mask, and are they disjoint?
    assert set(pixels_in_mask) == pixels_in_BASS_MzLS.union(pixels_in_DECaLS_NGC).union(pixels_in_DECaLS_SGC)
    assert pixels_in_BASS_MzLS.isdisjoint(pixels_in_DECaLS_NGC)
    assert pixels_in_BASS_MzLS.isdisjoint(pixels_in_DECaLS_SGC)
    assert pixels_in_DECaLS_NGC.isdisjoint(pixels_in_DECaLS_SGC)

    # Convert to np.array to finally save them into files
    pixels_in_BASS_MzLS = np.array(list(pixels_in_BASS_MzLS))
    pixels_in_DECaLS_NGC = np.array(list(pixels_in_DECaLS_NGC))
    pixels_in_DECaLS_SGC = np.array(list(pixels_in_DECaLS_SGC))

    # Sanity check: when plotting these masks together,
    # do they touch nicely at the borders?
    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_BASS_MzLS] = 1
    m[pixels_in_DECaLS_NGC] = 2
    m[pixels_in_DECaLS_SGC] = 3
    m[[32808]] = 4 # Quick trick: make one point near DEC=+90 have a larger value, to make colors in masks show as desired by cmap
    hp.mollview(m, nest=True, rot=[120, 0], title="BASS/MzLS and DECaLS masks together", cmap=cmap, cbar=False)
    hp.graticule()
    plt.show()

    # Print number of pixels in masks
    print("Number of pixels in BASS/MzLS mask", len(set(pixels_in_BASS_MzLS)))
    print("Number of pixels in DECaLS NGC mask", len(set(pixels_in_DECaLS_NGC)))
    print("Number of pixels in DECaLS SGC mask", len(set(pixels_in_DECaLS_SGC)))

    # Finally, save the resulting masks as bit arrays
    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_BASS_MzLS] = 1
    with open ("/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy", "wb") as f:
        np.save(f, m)

    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_DECaLS_NGC] = 1
    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_NGC_mask.npy", "wb") as f:
        np.save(f, m)

    m = np.zeros(NPIX, dtype=int)
    m[pixels_in_DECaLS_SGC] = 1
    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_SGC_mask.npy", "wb") as f:
        np.save(f, m)
