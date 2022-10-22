from astropy.io import fits
from desitarget.geomask import hp_in_box
import healpy as hp
import numpy as np
from tqdm import tqdm

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
    with fits.open(f"/cluster/scratch/lmachado/DESIImaging/dr9/{region}/footprint/{filename}") as f:
        data = f[1].data

        for brick in tqdm(data):
            radecbox = [brick["ra1"], brick["ra2"], brick["dec1"], brick["dec2"]]
            indices = hp_in_box(nside, radecbox, inclusive=True, fact=16)
            pixel_indices[region].update(indices)

    with open(f"/cluster/scratch/lmachado/DataProducts/footprint/pixels_in_{region}_footprint_for_nside_{nside}.npy", "wb") as f:
        np.save(f, np.array(list(pixel_indices[region])))
