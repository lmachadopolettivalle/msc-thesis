import numpy as np

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

from load_processed_target_data import load_processed_target_data
from mask import mask

REGION = "north"

# Load target data
print("Loading target data")
targets_nside = 512
targets = load_processed_target_data(region=REGION, extinction_correction=True)

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


# TODO remove me
# Use small number of objects to test code
targets = {
    k: v[:1000]
    for k, v in targets.items()
}
randoms = {
    k: v[:10000]
    for k, v in randoms.items()
}


# Apply mask to both targets and randoms
print("Applying mask to targets and randoms")
targets_masked = mask(targets["HPXPIXEL"], targets_nside, region=REGION)
targets_ids_in_mask = np.where(targets_masked > 0)[0]

targets = {
    k: v[targets_ids_in_mask]
    for k, v in targets.items()
}

randoms_masked = mask(randoms["HPXPIXEL"], randoms_nside, region=REGION)
randoms_ids_in_mask = np.where(randoms_masked > 0)[0]

randoms = {
    k: v[randoms_ids_in_mask]
    for k, v in randoms.items()
}

# Count number of targets and randoms
targets_count = len(targets["RA"])
randoms_count = len(randoms["RA"])

print(targets_count, randoms_count)

# Compute correlation function
print("Computing 2PCF")
nbins = 15
bins = np.logspace(np.log10(0.01), np.log10(10.0), nbins + 1)
nthreads = 1

DD_counts = DDtheta_mocks(1, nthreads, bins, targets["RA"], targets["DEC"])
RR_counts = DDtheta_mocks(1, nthreads, bins, randoms["RA"], randoms["DEC"])
DR_counts = DDtheta_mocks(0, nthreads, bins, targets["RA"], targets["DEC"], RA2=randoms["RA"], DEC2=randoms["DEC"])

wtheta = convert_3d_counts_to_cf(
    targets_count, targets_count,
    randoms_count, randoms_count,
    DD_counts, DR_counts,
    DR_counts, RR_counts
)

with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_bins.npy", "wb") as f:
    np.save(f, bins[:-1])
with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_wtheta.npy", "wb") as f:
    np.save(f, wtheta)
