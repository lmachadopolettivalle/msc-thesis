from astropy.io import fits
import numpy as np
import os
from tqdm import tqdm

REGION = "north"

def count_original_targets(region=REGION):
    print(f"Data for region: {region}")

    PATH_TO_SWEEP_FILES = f"/cluster/scratch/lmachado/DESIImaging/dr9/{region}/sweeps/"

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

def print_mask_areas():
    import healpy as hp

    NSIDE = 64
    NPIX = hp.nside2npix(NSIDE)
    PIXEL_AREA = hp.pixelfunc.nside2pixarea(NSIDE, degrees=True)
    SKY_AREA = NPIX * PIXEL_AREA

    with open ("/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy", "rb") as f:
        bass_mask = np.load(f)

    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_NGC_mask.npy", "rb") as f:
        decals_ngc_mask = np.load(f)

    with open ("/cluster/scratch/lmachado/DataProducts/masks/DECaLS_SGC_mask.npy", "rb") as f:
        decals_sgc_mask = np.load(f)

    bass_pixels = np.sum(bass_mask)
    decals_ngc_pixels = np.sum(decals_ngc_mask)
    decals_sgc_pixels = np.sum(decals_sgc_mask)

    print(f"Total number of pixels in sky at nside {NSIDE} is {NPIX}")
    print(f"Number of pixels in BASS/MzLS mask: {bass_pixels}")
    print(f"Number of pixels in DECaLS-NGC mask: {decals_ngc_pixels}")
    print(f"Number of pixels in DECaLS-SGC mask: {decals_sgc_pixels}")

    bass_area = bass_pixels * PIXEL_AREA
    decals_ngc_area = decals_ngc_pixels * PIXEL_AREA
    decals_sgc_area = decals_sgc_pixels * PIXEL_AREA

    print(f"BASS/MzLS mask area = {bass_area} deg^2, {bass_area / SKY_AREA} of the sky")
    print(f"DECaLS-NGC mask area = {decals_ngc_area} deg^2, {decals_ngc_area / SKY_AREA} of the sky")
    print(f"DECaLS-SGC mask area = {decals_sgc_area} deg^2, {decals_sgc_area / SKY_AREA} of the sky")

def print_target_morphtype_counts(region=REGION):
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_MORPHTYPE.npy", "rb") as f:
        morphtypes = np.load(f)

    print(f"Targets in region {region} BEFORE MASKING: {len(morphtypes)}")

    unique, unique_counts = np.unique(morphtypes, return_counts=True)
    print(f"Region {region}: {unique} | {unique_counts}")

if __name__ == "__main__":
    print_mask_areas()
    print_target_morphtype_counts(region="north")
    print_target_morphtype_counts(region="south")
