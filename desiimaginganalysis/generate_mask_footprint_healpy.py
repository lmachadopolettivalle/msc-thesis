from astropy.io import fits
from desitarget.geomask import hp_in_box
import numpy as np
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
