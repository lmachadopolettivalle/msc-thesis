from astropy.io import fits

from tqdm import tqdm

from desitarget.geomask import hp_in_box

import numpy as np

PATH_TO_FOOTPRINT_FILES = "/cluster/home/lmachado/msc-thesis/dr9/footprints/"
FOOTPRINT_FILES = ["survey-bricks-dr9-north.fits", "survey-bricks-dr9-south.fits"]

footprint = []

nside = 256

pixel_indices = set()

for filename in FOOTPRINT_FILES:
    with fits.open(f"{PATH_TO_FOOTPRINT_FILES}/{filename}") as f:
        data = f[1].data

        for brick in tqdm(data):
            radecbox = [brick["ra1"], brick["ra2"], brick["dec1"], brick["dec2"]]
            indices = hp_in_box(nside, radecbox, inclusive=True, fact=4)
            pixel_indices.update(indices)
    print(pixel_indices)

with open(f"pixels_in_footprint_for_nside_{nside}_method2.npy", "wb") as f:
    np.save(f, np.array(list(pixel_indices)))
