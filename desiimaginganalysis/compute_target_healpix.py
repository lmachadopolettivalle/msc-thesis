import healpy as hp
import numpy as np

REGION = "south"

NSIDE = 512

with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_RA.npy", "rb") as f:
    RAs = np.load(f)
with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_DEC.npy", "rb") as f:
    DECs = np.load(f)

pixel_ids = hp.pixelfunc.ang2pix(
    NSIDE,
    RAs, # 0 to 360
    DECs, # -90 to +90
    nest=True,
    lonlat=True,
)

with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_HPXPIXEL_HPXNSIDE_{NSIDE}.npy", "wb") as f:
    np.save(f, np.array(pixel_ids))
