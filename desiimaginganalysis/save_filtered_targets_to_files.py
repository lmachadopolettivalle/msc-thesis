import healpy as hp
import numpy as np

from select_imaging_targets import targets, REGION
print("Done selecting targets.")
print("About to start filtering columns and saving to files.")

# Desired columns to be stored
COLUMNS = [
    "RELEASE", "BRICKID", "BRICKNAME", "BRICK_OBJID", "MORPHTYPE",
    "RA", "RA_IVAR", "DEC", "DEC_IVAR",
    "FLUX_G", "FLUX_R", "FLUX_Z",
    "FLUX_IVAR_G", "FLUX_IVAR_R", "FLUX_IVAR_Z",
    "MW_TRANSMISSION_G", "MW_TRANSMISSION_R", "MW_TRANSMISSION_Z",
    "PHOTSYS", "TARGETID",
    "DESI_TARGET", "BGS_TARGET",
]

data = targets[:][COLUMNS]

# Save one column per file to use memory more efficiently

for column in COLUMNS:
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_{column}.npy", "wb") as f:
        np.save(f, data[:][column])

# Finally, compute target HEALPix pixel IDs

NSIDE = 512

pixel_ids = hp.pixelfunc.ang2pix(
    NSIDE,
    data[:]["RA"], # 0 to 360
    data[:]["DEC"], # -90 to +90
    nest=True,
    lonlat=True,
)

with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_HPXPIXEL_HPXNSIDE_{NSIDE}.npy", "wb") as f:
    np.save(f, np.array(pixel_ids))
