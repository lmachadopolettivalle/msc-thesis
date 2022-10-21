import healpy as hp
from matplotlib import pyplot as plt
import numpy as np

from desitarget.io import desitarget_resolve_dec
from desitarget.geomask import nside2nside

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

# TODO uncomment south targets
south_target_pixels = south_footprint_pixels.copy()
#with open(SOUTH_TARGET_PIXELS_FILE, "rb") as f:
#    south_target_pixels = np.load(f)



north_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, north_footprint_pixels)
north_target_pixels = nside2nside(NSIDE_FILES, NSIDE, north_target_pixels)
south_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, south_footprint_pixels)
south_target_pixels = nside2nside(NSIDE_FILES, NSIDE, south_target_pixels)

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
hp.mollview(m, nest=True, rot=[120, 0])
hp.graticule()
plt.show()


############
# North mask
############
print("NORTH")

# Get all neighbor pixels of all target pixels
# Note: since phi=None, theta is taken as a list of pixel IDs
all_neighbor_pixels = hp.pixelfunc.get_all_neighbours(
    NSIDE,
    theta=north_target_pixels,
    phi=None,
    nest=True,
)

all_neighbor_pixels = all_neighbor_pixels.flatten()

# Remove target pixels from the neighbor list
right_outside_pixels = np.array(list(
    set(all_neighbor_pixels).difference(set(north_target_pixels))
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
pixels_in_north_mask = set(north_target_pixels).difference(set(right_outside_neighbors))

# Take intersection with northern region to remove any boundary pixels
pixels_in_north_mask = np.array(list(
    pixels_in_north_mask.intersection(set(nest_north_pixels))
))
pixels_in_north_mask = np.array(list(pixels_in_north_mask))

# Visualize mask and targets together at NSIDE
m = np.zeros(NPIX)
m[north_target_pixels] = 1
m[pixels_in_north_mask] = 2

hp.mollview(m, nest=True, rot=[120, 0])
plt.show()


# Sanity check:
# Is the mask totally inside the footprint?
mask_footprint_intersection = set(pixels_in_north_mask).intersection(set(north_footprint_pixels))
print("Intersection between mask and footprint:", len(mask_footprint_intersection))
print("Pixels in mask:", len(set(pixels_in_north_mask)))
# If the lenghts are not the same, the mask is not entirely contained in the footprint
assert len(mask_footprint_intersection) == len(set(pixels_in_north_mask))

# Sanity check:
# Is the mask totally inside the targets area?
mask_target_intersection = set(pixels_in_north_mask).intersection(set(north_target_pixels))
print("Intersection between mask and targets:", len(mask_target_intersection))
print("Pixels in mask:", len(set(pixels_in_north_mask)))
m = np.zeros(NPIX)
m[np.array(list(set(pixels_in_north_mask).difference(set(north_target_pixels))), dtype=int)] = 1
hp.mollview(m, nest=True, rot=[120, 0])
hp.graticule()
plt.show()

# Sanity check:
# Are there any holes in the gray footprint?
footprint_holes = set(pixels_in_north_mask).difference(set(north_footprint_pixels))
assert len(footprint_holes) == 0

# Sanity check:
# Are there any holes in the targets pixels?
target_holes = set(pixels_in_north_mask).difference(set(north_target_pixels))
print("Holes in target pixels (i.e. number of pixels in mask but not in target):", len(target_holes))

m = np.zeros(NPIX)
m[pixels_in_north_mask] = 1
hp.mollview(m, nest=True, rot=[120, 0])
hp.graticule()
plt.show()


# Compute area percentages
footprint_loss = (len(set(north_footprint_pixels)) - len(set(pixels_in_north_mask))) / len(set(north_footprint_pixels))
target_loss = (len(set(north_target_pixels)) - len(set(pixels_in_north_mask))) / len(set(north_target_pixels))

print("Gray footprint loss:", footprint_loss)
print("Target area loss:", target_loss)


############
# South mask
############
print("SOUTH")

# Get all neighbor pixels of all target pixels
# Note: since phi=None, theta is taken as a list of pixel IDs
all_neighbor_pixels = hp.pixelfunc.get_all_neighbours(
    NSIDE,
    theta=south_target_pixels,
    phi=None,
    nest=True,
)

all_neighbor_pixels = all_neighbor_pixels.flatten()

# Remove target pixels from the neighbor list
right_outside_pixels = np.array(list(
    set(all_neighbor_pixels).difference(set(south_target_pixels))
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
pixels_in_south_mask = set(south_target_pixels).difference(set(right_outside_neighbors))

# Take intersection with southern region to remove any boundary pixels
print("HAHAHA")
print(len(pixels_in_south_mask))
pixels_in_south_mask = np.array(list(
    pixels_in_south_mask.intersection(set(nest_south_pixels))
))
print(len(pixels_in_south_mask))
print("HAHAHA")

pixels_in_south_mask = np.array(list(pixels_in_south_mask))

# Visualize mask and targets together at NSIDE
m = np.zeros(NPIX)
m[south_target_pixels] = 1
m[pixels_in_south_mask] = 2

hp.mollview(m, nest=True, rot=[120, 0])
plt.show()



# Sanity check:
# Is the mask totally inside the footprint?
mask_footprint_intersection = set(pixels_in_south_mask).intersection(set(south_footprint_pixels))
print("Intersection between mask and footprint:", len(mask_footprint_intersection))
print("Pixels in mask:", len(set(pixels_in_south_mask)))
# If the lenghts are not the same, the mask is not entirely contained in the footprint
assert len(mask_footprint_intersection) == len(set(pixels_in_south_mask))

"""
# TODO include this code when we have south targets
# Sanity check:
# Is the mask totally inside the targets area?
mask_target_intersection = set(pixels_in_south_mask).intersection(set(south_target_pixels))
print("Intersection between mask and targets:", len(mask_target_intersection))
print("Pixels in mask:", len(set(pixels_in_south_mask)))
m = np.zeros(NPIX)
m[np.array(list(set(pixels_in_south_mask).difference(set(south_target_pixels))))] = 1
hp.mollview(m, nest=True, rot=[120, 0])
hp.graticule()
plt.show()
"""

# Sanity check:
# Are there any holes in the gray footprint?
footprint_holes = set(pixels_in_south_mask).difference(set(south_footprint_pixels))
assert len(footprint_holes) == 0

# Sanity check:
# Are there any holes in the targets pixels?
target_holes = set(pixels_in_south_mask).difference(set(south_target_pixels))
print("Holes in target pixels (i.e. number of pixels in mask but not in target):", len(target_holes))

m = np.zeros(NPIX)
m[pixels_in_south_mask] = 1
hp.mollview(m, nest=True, rot=[120, 0])
hp.graticule()
plt.show()


# Compute area percentages
footprint_loss = (len(set(south_footprint_pixels)) - len(set(pixels_in_south_mask))) / len(set(south_footprint_pixels))
#target_loss = (len(set(south_target_pixels)) - len(set(pixels_in_south_mask))) / len(set(south_target_pixels))

print("Gray footprint loss:", footprint_loss)
#print("Target area loss", target_loss)

