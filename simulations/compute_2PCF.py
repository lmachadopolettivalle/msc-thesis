# Using x, y, z coordinates from galaxies
# as outputted by SHAM,
# determine their RA and Dec,
# and compute the angular correlation function

import argparse

import numpy as np
import os

from desitarget.io import desitarget_resolve_dec

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

import sys
sys.path.append("..")
from desiimaginganalysis.mask import mask

import directories

from load_sham_galaxies import load_sham_galaxies

from sham_model_constants import BLUE, RED

# Equation 2 in Zarrouk+21
def rmag_bass_correction(unprimed_rmags, unprimed_gmags):
    assert len(unprimed_rmags) == len(unprimed_gmags)
    return unprimed_rmags - 0.039 * (unprimed_gmags - unprimed_rmags) + 0.011

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord
BANDS = ["mag_g", "mag_r", "mag_z"]

# Whether to use r-primed
# Only applies to BASS_MzLS objects
USE_MAG_R = True

parser = argparse.ArgumentParser()
parser.add_argument("--run_id", type=int, required=True)
parser.add_argument("--region", type=str, required=True)
DEFAULT_SEED = "0"
parser.add_argument("--seed", type=str, required=False, default=DEFAULT_SEED)

args = parser.parse_args()
run_id = args.run_id
DESI_region = args.region
seed = args.seed
if seed == DEFAULT_SEED:
    seed = None

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

# Path where to save 2PCF computed values
PATH_2PCF = directories.path_2PCF(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_region,
    run_id=run_id,
    seed=seed,
)
if os.path.isdir(PATH_2PCF):
    print(f"{PATH_2PCF} directory already exists.")
else:
    print(f"Creating new output directory, {PATH_2PCF} ...")
    os.mkdir(PATH_2PCF)
    print("Created output directory successfully.")

print("Loading galaxy data")
# Load galaxy data from SHAM
galaxies = load_sham_galaxies(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_region,
    run_id=run_id,
    seed=seed,
)

# For BASS, compute r-primed magnitudes
if (USE_MAG_R is True):
    if (DESI_region == directories.BASS_MzLS):
        galaxies["mag_r"] = rmag_bass_correction(
            galaxies["mag_r"],
            galaxies["mag_g"],
        )
    elif (DESI_region == directories.FULLSKY):
        dec_mask = (galaxies["DEC"] >= desitarget_resolve_dec())
        galaxies["mag_r"][dec_mask] = rmag_bass_correction(
            galaxies["mag_r"][dec_mask],
            galaxies["mag_g"][dec_mask],
        )

print("Done loading galaxy data")

# Load randoms, apply mask
print("Loading randoms")
randoms_nside = 256
randoms = {}
with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_RA.npy", "rb") as f:
    randoms["RA"] = np.load(f)
with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_DEC.npy", "rb") as f:
    randoms["DEC"] = np.load(f)
with open(f"/cluster/scratch/lmachado/DataProducts/randoms/randoms_pixels_NSIDE_{randoms_nside}.npy", "rb") as f:
    randoms["HPXPIXEL"]= np.load(f)

# Range used to filter randoms was chosen via trial and error, to make sure
# there is a similar number of randoms and targets
RANDOMS_SIZE = 10000000
if DESI_region == directories.FULLSKY:
    RANDOMS_SIZE = 20000000
randoms = {
    k: v[:RANDOMS_SIZE]
    for k, v in randoms.items()
}

# Apply mask to randoms
print("Applying mask to randoms")

if DESI_region == directories.FULLSKY:
    REGIONS = (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC, )
else:
    REGIONS = (DESI_region, )


randoms_masked = mask(randoms["HPXPIXEL"], randoms_nside, regions=REGIONS)
randoms_ids_in_mask = np.where(randoms_masked > 0)[0]

randoms = {
    k: v[randoms_ids_in_mask]
    for k, v in randoms.items()
}

print("Total random count (after masking):", len(randoms["RA"]))

# Compute 2PCF for different r-mag bins
rmag_bins = [
    [14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

# Compute correlation function

# Many bins, useful for plotting
nbins = 24
bins = np.logspace(np.log10(1e-3), np.log10(20), nbins + 1, endpoint=True)

"""
# Few bins, useful for covariance matrix
nbins = 7
tmp_bins = np.logspace(np.log10(6e-2), np.log10(3), nbins)
bins = (tmp_bins[:-1] + tmp_bins[1:]) / 2
"""


nthreads = 8

# Pre-compute randoms-related data
print("Computing mocks for randoms")
RR_counts = DDtheta_mocks(1, nthreads, bins, randoms["RA"], randoms["DEC"])
randoms_count = len(randoms["RA"])

print("Computing 2PCF for regions", REGIONS)


# Compute 2PCF, for blue and red galaxies separately
for color_name, color_value in (
        ("total", None),
        #("blue", BLUE),
        #("red", RED)
    ):
    print(f"Computing 2PCF for {color_name} galaxies")

    if color_value is None:
        # For total, use all galaxies when computing 2PCF
        color_galaxies = galaxies
    else:
        color_bitmask = np.where(galaxies["blue_red"] == color_value)[0]
        color_galaxies = {
            k: v[color_bitmask]
            for k, v in galaxies.items()
        }


    for rmag_low, rmag_high in rmag_bins:
        print(rmag_low, rmag_high)
        rmag_ids = np.where(
            (rmag_low <= color_galaxies["mag_r"]) &
            (color_galaxies["mag_r"] < rmag_high)
        )[0]
        rmag_filtered_targets = {
            k: v[rmag_ids]
            for k, v in color_galaxies.items()
        }

        # Count number of targets and randoms
        targets_count = len(rmag_filtered_targets["RA"])

        print("Targets count:", targets_count)

        DD_counts = DDtheta_mocks(1, nthreads, bins, rmag_filtered_targets["RA"], rmag_filtered_targets["DEC"])
        DR_counts = DDtheta_mocks(0, nthreads, bins, rmag_filtered_targets["RA"], rmag_filtered_targets["DEC"], RA2=randoms["RA"], DEC2=randoms["DEC"])

        wtheta = convert_3d_counts_to_cf(
            targets_count, targets_count,
            randoms_count, randoms_count,
            DD_counts, DR_counts,
            DR_counts, RR_counts
        )

        with open(f"{PATH_2PCF}/simulated_{color_name}_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy", "wb") as f:
            np.save(f, bins[:-1])
        with open(f"{PATH_2PCF}/simulated_{color_name}_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy", "wb") as f:
            np.save(f, wtheta)

print("Done computing 2PCF")
