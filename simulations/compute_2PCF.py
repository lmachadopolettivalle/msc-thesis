# Using x, y, z coordinates from galaxies
# as outputted by SHAM,
# determine their RA and Dec,
# and compute the angular correlation function

import numpy as np

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

import sys
sys.path.append("..")
from desiimaginganalysis.mask import mask

from sham_model_constants import BLUE, RED

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord.
BANDS = ["mag_g", "mag_r", "mag_z"]

BASS_MzLS = "BASS-MzLS"
REGIONS = (BASS_MzLS, )

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 512

# Load x, y, z positions
# Path to output data from SHAM
SHAM_OUTPUT_PATH = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{PARTICLE_COUNT_PINOCCHIO}cubed/interpolation_outputs/"

galaxies = {}
for coord in ("x_coord", "y_coord", "z_coord"):
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_{coord}.npy"

    galaxies[coord] = np.load(filename)

# Load magnitudes
for band in BANDS:
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_app_{band}.npy"

    galaxies[band] = np.load(filename)


# Load whether galaxies are blue or red
filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_blue_red.npy"

galaxies["blue_red"] = np.load(filename)

# Convert 3D positions into RA, Dec
radii = np.sqrt(galaxies["x_coord"]**2 + galaxies["y_coord"]**2 + galaxies["z_coord"]**2)

theta = np.arccos(galaxies["z_coord"] / radii)
phi = np.arctan2(galaxies["y_coord"], galaxies["x_coord"])

# Note that phi is in range [-pi, pi], but for healpy, must be in range [0, 360 degrees]
phi[phi < 0] += 2 * np.pi

galaxies["RA"] = np.degrees(phi)
galaxies["DEC"] = np.degrees(np.pi/2 - theta)

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
randoms = {
    k: v[:3000000]
    for k, v in randoms.items()
}

# Apply mask to randoms
print("Applying mask to randoms")

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
nbins = 24
bins = np.logspace(np.log10(1e-3), np.log10(20), nbins + 1, endpoint=True)
nthreads = 8

# Pre-compute randoms-related data
print("Computing mocks for randoms")
RR_counts = DDtheta_mocks(1, nthreads, bins, randoms["RA"], randoms["DEC"])
randoms_count = len(randoms["RA"])

print("Computing 2PCF for regions", REGIONS)


# Compute 2PCF, for blue and red galaxies separately
for color_name, color_value in (("blue", BLUE), ("red", RED)):
    print(f"Computing 2PCF for color {color_name}")
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

        with open(f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy", "wb") as f:
            np.save(f, bins[:-1])
        with open(f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy", "wb") as f:
            np.save(f, wtheta)

print("Done computing 2PCF")
