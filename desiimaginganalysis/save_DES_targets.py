import time

from load_processed_target_data import load_processed_target_data
from constants import DES

# Visualize DES mask
DESMASK = hp.fitsfunc.read_map("/cluster/work/refregier/bernerp/DES/DES_Y1/y1a1_gold_1.0.2_wide_footprint_4096.fit", nest=False)
hp.mollview(DESMASK, nest=False, rot=[120, 0])
plt.show()

# Load target data with DES mask
start = time.perf_counter()
DES_targets = load_processed_target_data(regions=[DES], extinction_correction=True, apply_mask=True)
end = time.perf_counter()

print(end-start)

# Compute raw magnitudes
DES_targets["MAG_R_RAW"] = 22.5 - 2.5*np.log10(DES_targets["FLUX_R"].clip(1e-16))
DES_targets["MAG_G_RAW"] = 22.5 - 2.5*np.log10(DES_targets["FLUX_G"].clip(1e-16))
DES_targets["MAG_Z_RAW"] = 22.5 - 2.5*np.log10(DES_targets["FLUX_Z"].clip(1e-16))

# Make histogram of difference in r magnitudes
plt.hist(DES_targets["MAG_R"] - DES_targets["MAG_R_RAW"], bins=100); plt.xlabel("r (corrected) - r (raw)"); plt.ylabel("Count"); plt.show()

# Save DES-selected catalogs
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_RA.npy", "wb") as f:
    np.save(f, DES_targets["RA"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_DEC.npy", "wb") as f:
    np.save(f, DES_targets["DEC"])

with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_G_RAW.npy", "wb") as f:
    np.save(f, DES_targets["MAG_G_RAW"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_R_RAW.npy", "wb") as f:
    np.save(f, DES_targets["MAG_R_RAW"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_Z_RAW.npy", "wb") as f:
    np.save(f, DES_targets["MAG_Z_RAW"])

with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_G_CORRECTED.npy", "wb") as f:
    np.save(f, DES_targets["MAG_G"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_R_CORRECTED.npy", "wb") as f:
    np.save(f, DES_targets["MAG_R"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MAG_Z_CORRECTED.npy", "wb") as f:
    np.save(f, DES_targets["MAG_Z"])

with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MW_TRANSMISSION_G.npy", "wb") as f:
    np.save(f, DES_targets["MW_TRANSMISSION_G"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MW_TRANSMISSION_R.npy", "wb") as f:
    np.save(f, DES_targets["MW_TRANSMISSION_R"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_MW_TRANSMISSION_Z.npy", "wb") as f:
    np.save(f, DES_targets["MW_TRANSMISSION_Z"])

with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_FLUX_G.npy", "wb") as f:
    np.save(f, DES_targets["FLUX_G"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_FLUX_R.npy", "wb") as f:
    np.save(f, DES_targets["FLUX_R"])
with open("/cluster/scratch/lmachado/DataProducts/targets/DES/targets_FLUX_Z.npy", "wb") as f:
    np.save(f, DES_targets["FLUX_Z"])
