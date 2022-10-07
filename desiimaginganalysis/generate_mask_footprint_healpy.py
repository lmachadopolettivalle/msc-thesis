from astropy.io import fits
from desitarget.geomask import hp_in_box
import healpy as hp
import numpy as np

def mask(ras, decs, nside=512, nest=True, region="both"):
    # region can be "north", "south", or "both"
    # ras must be a list of RA values between 0 and 360
    # decs must be a list of Dec values between -90 and +90
    # ras and decs must have same length
    # returns: list of booleans, indicating whether each target is in the footprint

    # Example call:
    # ras = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 210, 240, 270]
    # decs = [15]*len(ras)
    # m = mask(ras, decs, nside=512, region='both')
    # print(m)
    #   array([ True,  True,  True, False, False, False, False, False,  True,
    #   True,  True,  True,  True,  True,  True,  True])

    # Load footprints
    if (region == "north") or (region == "both"):
        with open(f"pixels_in_north_footprint_for_nside_{nside}.npy", "rb") as f:
            footprint_north = np.load(f)
    if (region == "south") or (region == "both"):
        with open(f"pixels_in_south_footprint_for_nside_{nside}.npy", "rb") as f:
            footprint_south = np.load(f)

    footprint = set()
    if (region == "north") or (region == "both"):
        footprint.update(footprint_north)
    if (region == "south") or (region == "both"):
        footprint.update(footprint_south)

    footprint = np.array(list(footprint))

    pixel_ids = hp.ang2pix(nside, ras, decs, nest=nest, lonlat=True)

    return np.isin(pixel_ids, footprint)


if __name__ == "__main__":
    from tqdm import tqdm

    PATH_TO_FOOTPRINT_FILES = "/cluster/home/lmachado/msc-thesis/dr9/footprints/"
    FOOTPRINT_FILES = {
        "north": "survey-bricks-dr9-north.fits",
        "south": "survey-bricks-dr9-south.fits",
    }

    pixel_indices = {
        "north": set(),
        "south": set(),
    }

    nside = 512

    for region, filename in FOOTPRINT_FILES.items():
        with fits.open(f"{PATH_TO_FOOTPRINT_FILES}/{filename}") as f:
            data = f[1].data

            for brick in tqdm(data):
                radecbox = [brick["ra1"], brick["ra2"], brick["dec1"], brick["dec2"]]
                indices = hp_in_box(nside, radecbox, inclusive=True, fact=16)
                pixel_indices[region].update(indices)

        with open(f"pixels_in_{region}_footprint_for_nside_{nside}.npy", "wb") as f:
            np.save(f, np.array(list(pixel_indices[region])))
