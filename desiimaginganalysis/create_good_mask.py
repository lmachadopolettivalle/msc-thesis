import healpy as hp
from matplotlib import pyplot as plt
import numpy as np

# Load footprint and targets
NSIDE_FILES = 512
NSIDE = 256

NPIX = hp.nside2npix(NSIDE)

NORTH_FOOTPRINT_FILE = f"./DataProducts/footprint/pixels_in_north_footprint_for_nside_{NSIDE_FILES}.npy"
SOUTH_FOOTPRINT_FILE = f"./DataProducts/footprint/pixels_in_south_footprint_for_nside_{NSIDE_FILES}.npy"

NORTH_TARGET_PIXELS_FILE = f"./DataProducts/targets/north/targets_HPXPIXEL_HPXNSIDE_{NSIDE_FILES}.npy"
SOUTH_TARGET_PIXELS_FILE = f"./DataProducts/targets/south/targets_HPXPIXEL_HPXNSIDE_{NSIDE_FILES}.npy"

with open(NORTH_FOOTPRINT_FILE, "rb") as f:
    north_footprint_pixels = np.load(f)
with open(NORTH_TARGET_PIXELS_FILE, "rb") as f:
    north_target_pixels = np.load(f)
with open(SOUTH_FOOTPRINT_FILE, "rb") as f:
    south_footprint_pixels = np.load(f)
# TODO uncomment south targets
south_target_pixels = np.array([])
#with open(SOUTH_TARGET_PIXELS_FILE, "rb") as f:
#    south_target_pixels = np.load(f)

# Test converting to smaller nside for gaps in target map
def nside2nside(nside, nsidenew, pixlist):
    """Change a list of HEALPixel numbers to a different NSIDE.

    Parameters
    ----------
    nside : :class:`int`
        The current HEALPixel nside number (NESTED scheme).
    nsidenew : :class:`int`
        The new HEALPixel nside number (NESTED scheme).
    pixlist : :class:`list` or `~numpy.ndarray`
        The list of HEALPixels to be changed.

    Returns
    -------
    :class:`~numpy.ndarray`
        The altered list of HEALPixels.

    Notes
    -----
        - The size of the input list will be altered. For instance,
          nside=2, pixlist=[0,1] is covered by only pixel [0] at
          nside=1 but by pixels [0, 1, 2, 3, 4, 5, 6, 7] at nside=4.
        - Doesn't check that the passed pixels are valid at `nside`.
    """
    # ADM sanity check that we're in the nested scheme.
    #check_nside([nside, nsidenew])

    pixlist = np.atleast_1d(pixlist)

    # ADM remember to use integer division throughout.
    # ADM if nsidenew is smaller (at a lower resolution), then
    # ADM downgrade the passed pixel numbers.
    if nsidenew <= nside:
        fac = (nside//nsidenew)**2
        pixlistnew = np.array(list(set(pixlist//fac)))
    else:
        # ADM if nsidenew is larger (at a higher resolution), then
        # ADM upgrade the passed pixel numbers.
        fac = (nsidenew//nside)**2
        pixlistnew = []
        for pix in pixlist:
            pixlistnew.append(np.arange(pix*fac, pix*fac+fac))
        pixlistnew = np.hstack(pixlistnew)

    return pixlistnew

north_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, north_footprint_pixels)
north_target_pixels = nside2nside(NSIDE_FILES, NSIDE, north_target_pixels)
south_footprint_pixels = nside2nside(NSIDE_FILES, NSIDE, south_footprint_pixels)
south_target_pixels = nside2nside(NSIDE_FILES, NSIDE, south_target_pixels)



# Obtain pixels in correct region

# TODO replace with 
#from desitarget.io import desitarget_resolve_dec
def desitarget_resolve_dec():
    return 32.375

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
"""
m = np.zeros(NPIX)
m[north_footprint_pixels] = 1
m[south_footprint_pixels] = 2
hp.mollview(m, nest=True, rot=[180, 0])
hp.graticule()
plt.show()
"""


# Focusing on North,
# find a polygon that covers a large enough part of the footprint

# TODO remove this
def hp_in_dec_range(nside, decmin, decmax, inclusive=True):
    """HEALPixels in a specified range of Declination.

    Parameters
    ----------
    nside : :class:`int`
        (NESTED) HEALPixel nside.
    decmin, decmax : :class:`float`
        Declination range (degrees).
    inclusive : :class:`bool`, optional, defaults to ``True``
        see documentation for `healpy.query_strip()`.

    Returns
    -------
    :class:`list`
        (Nested) HEALPixels at `nside` in the specified Dec range.

    Notes
    -----
        - Just syntactic sugar around `healpy.query_strip()`.
        - `healpy.query_strip()` isn't implemented for the NESTED scheme
          in early healpy versions, so this queries in the RING scheme
          and then converts to the NESTED scheme.
    """
    # ADM convert Dec to co-latitude in radians.
    # ADM remember that, min/max swap because of the -ve sign.
    thetamin = np.radians(90.-decmax)
    thetamax = np.radians(90.-decmin)

    # ADM determine the pixels that touch the box.
    pixring = hp.query_strip(nside, thetamin, thetamax,
                             inclusive=inclusive, nest=False)
    pixnest = hp.ring2nest(nside, pixring)

    return pixnest

# TODO replace with from desitarget.geomask import hp_in_box
def hp_in_box(nside, radecbox, inclusive=True, fact=4):
    """Determine which HEALPixels touch an RA, Dec box.

    Parameters
    ----------
    nside : :class:`int`
        (NESTED) HEALPixel nside.
    radecbox : :class:`list`
        4-entry list of coordinates [ramin, ramax, decmin, decmax]
        forming the edges of a box in RA/Dec (degrees).
    inclusive : :class:`bool`, optional, defaults to ``True``
        see documentation for `healpy.query_polygon()`.
    fact : :class:`int`, optional defaults to 4
        see documentation for `healpy.query_polygon()`.

    Returns
    -------
    :class:`list`
        HEALPixels at the passed `nside` that touch the RA/Dec box.

    Notes
    -----
        - Uses `healpy.query_polygon()` to retrieve the RA geodesics
          and then :func:`hp_in_dec_range()` to limit by Dec.
        - When the RA range exceeds 180o, `healpy.query_polygon()`
          defines the range as that with the smallest area (i.e the box
          can wrap-around in RA). To avoid any ambiguity, this function
          will only limit by the passed Decs in such cases.
        - Only strictly correct for Decs from -90+1e-3(o) to 90-1e3(o).
    """
    ramin, ramax, decmin, decmax = radecbox

    # ADM area enclosed isn't well-defined if RA covers more than 180o.
    if np.abs(ramax-ramin) <= 180.:
        # ADM retrieve RA range. The 1e-3 prevents edge effects near poles.
        npole, spole = 90-1e-3, -90+1e-3
        # ADM convert RA/Dec to co-latitude and longitude in radians.
        rapairs = np.array([ramin, ramin, ramax, ramax])
        decpairs = np.array([spole, npole, npole, spole])
        thetapairs, phipairs = np.radians(90.-decpairs), np.radians(rapairs)

        # ADM convert to Cartesian vectors remembering to transpose
        # ADM to pass the array to query_polygon in the correct order.
        vecs = hp.dir2vec(thetapairs, phipairs).T

        # ADM determine the pixels that touch the RA range.
        pixra = hp.query_polygon(nside, vecs,
                                 inclusive=inclusive, fact=fact, nest=True)
    else:
        log.warning('Max RA ({}) and Min RA ({}) separated by > 180o...'
                    .format(ramax, ramin))
        log.warning('...will only limit to passed Declinations'
                    .format(nside))
        pixra = np.arange(hp.nside2npix(nside))

    # ADM determine the pixels that touch the Dec range.
    pixdec = hp_in_dec_range(nside, decmin, decmax, inclusive=inclusive)

    # ADM return the pixels in the box.
    pixnum = list(set(pixra).intersection(set(pixdec)))

    return pixnum

pixels_in_north_mask = hp_in_box(
    NSIDE,
    radecbox=[105, 280, desitarget_resolve_dec()-2, 80],
    inclusive=False,
)

# Take intersection with northern region to remove any boundary pixels
pixels_in_north_mask = np.array(list(
    set(pixels_in_north_mask).intersection(set(nest_north_pixels))
))

# Sanity check:
# Is the mask totally inside the target footprint?
mask_footprint_intersection = set(pixels_in_north_mask).intersection(set(north_footprint_pixels))
print("Intersection between mask and footprint:", len(mask_footprint_intersection))
print("Pixels in mask:", len(set(pixels_in_north_mask)))
# If the lenghts are not the same, the mask is not entirely contained in the footprint
assert len(mask_footprint_intersection) == len(set(pixels_in_north_mask))

mask_target_intersection = set(pixels_in_north_mask).intersection(set(north_target_pixels))
print("Intersection between mask and targets:", len(mask_target_intersection))
print("Pixels in mask:", len(set(pixels_in_north_mask)))
m = np.zeros(NPIX)
m[np.array(list(set(pixels_in_north_mask).difference(set(north_target_pixels))))] = 1
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
print("Target area loss", target_loss)
