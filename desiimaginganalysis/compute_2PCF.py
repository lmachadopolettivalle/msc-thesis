import numpy as np

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

from load_processed_target_data import load_processed_target_data, BGS_BRIGHT, BGS_FAINT
from mask import mask

REGION = "north"

# Load target data
print("Loading target data")
targets_nside = 512
targets = load_processed_target_data(region=REGION, extinction_correction=True, apply_mask=True)

# Load randoms
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
    k: v[:40000000]
    for k, v in randoms.items()
}


# Apply mask to randoms
print("Applying mask to randoms")

randoms_masked = mask(randoms["HPXPIXEL"], randoms_nside, region=REGION)
randoms_ids_in_mask = np.where(randoms_masked > 0)[0]

randoms = {
    k: v[randoms_ids_in_mask]
    for k, v in randoms.items()
}

print("Total target count (after masking):", len(targets["RA"]))
print("Total random count (after masking):", len(randoms["RA"]))


# Compute 2PCF for each group of targets
bgs_category_target_ids = {
    "Bright": np.where(targets["BGS_TARGET"] == BGS_BRIGHT)[0],
    "Faint": np.where(targets["BGS_TARGET"] == BGS_FAINT)[0],
}

rmag_bins = {
    "Bright": [
        [15, 16],
        [16, 17],
        [17, 18],
        [18, 19],
        [19, 19.5],
    ],
    "Faint": [
        [19.5, 20],
        [20, 21],
    ],
}

# Compute correlation function
print("Computing 2PCF")
nbins = 25
bins = np.logspace(np.log10(1e-3), np.log10(10), nbins + 1)
nthreads = 2

# Pre-compute randoms-related data
RR_counts = DDtheta_mocks(1, nthreads, bins, randoms["RA"], randoms["DEC"])
randoms_count = len(randoms["RA"])

for bgs_category_name, ids in bgs_category_target_ids.items():
    print(bgs_category_name)
    bgs_category_filtered_targets = {
        k: v[ids]
        for k, v in targets.items()
    }

    for rmag_low, rmag_high in rmag_bins[bgs_category_name]:
        print(rmag_low, rmag_high)
        rmag_ids = np.where(
            (rmag_low <= bgs_category_filtered_targets["MAG_R"]) &
            (bgs_category_filtered_targets["MAG_R"] < rmag_high)
        )[0]
        rmag_filtered_targets = {
            k: v[rmag_ids]
            for k, v in bgs_category_filtered_targets.items()
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

        with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_bins_{bgs_category_name}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}.npy", "wb") as f:
            np.save(f, bins[:-1])
        with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_wtheta_{bgs_category_name}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}.npy", "wb") as f:
            np.save(f, wtheta)


print("Done computing 2PCF")
