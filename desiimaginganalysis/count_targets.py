from astropy.io import fits
import numpy as np
import os
from tqdm import tqdm

REGION = "north"
print(f"Data for region: {REGION}")

PATH_TO_SWEEP_FILES = f"/cluster/scratch/lmachado/DESIImaging/dr9/{REGION}/sweeps/"

FILENAMES = os.listdir(PATH_TO_SWEEP_FILES)

FILES = [f"{PATH_TO_SWEEP_FILES}/{filename}" for filename in FILENAMES]

# Galaxy types
TYPE_COUNTS = {
    k: 0
    for k in ("DEV", "DUP", "EXP", "PSF", "REX", "SER")
}

for filepath in tqdm(FILES):
    with fits.open(filepath) as f:
        data = f[1].data

        for k in TYPE_COUNTS.keys():
            TYPE_COUNTS[k] += np.sum(data["TYPE"] == k)

print(TYPE_COUNTS)
