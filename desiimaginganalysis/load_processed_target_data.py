# Read file with data from selected targets
import numpy as np
from desitarget.targets import bgs_mask

# Parameters for magnitude computation
CLIP_FLUX = 1e-16

BGS_BRIGHT = 1
BGS_FAINT = 0
# Note some targets may belong to WISE. We will not use these here.
BGS_WISE = -1

BANDS = ["g", "r", "z"]

def mag_from_flux(fluxes, mw_transmissions, extinction_correction=True):
    if extinction_correction:
        mags = 22.5 - 2.5 * np.log10((fluxes / mw_transmissions).clip(CLIP_FLUX))
    else:
        mags = 22.5 - 2.5 * np.log10(fluxes.clip(CLIP_FLUX))

    return mags

def remove_spurious_objects(target_dict):
    # Modify mags_dict to remove entries from all its mags,
    # if mag in at least one of them is due to clipping of the flux.
    # Explanation: due to the clip(1e-16), there will be a few objects (< 10 or so)
    # with magnitudes == 62.5.
    # These objects can be ignored.

    clip_mag = 22.5 - 2.5 * np.log10(CLIP_FLUX)

    r_idx = target_dict["MAG_R"] < clip_mag
    g_idx = target_dict["MAG_G"] < clip_mag
    z_idx = target_dict["MAG_Z"] < clip_mag

    all_idx = r_idx & g_idx & z_idx

    return {
        k: v[all_idx]
        for k, v in target_dict.items()
    }

def load_processed_target_data(region="north", extinction_correction=True):
    # Return dict, with values being arrays of the following fields for targets:
    # RA, DEC, Magnitudes, BGS Bright or Faint
    NSIDE = 512

    data = {}

    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_RA.npy", "rb") as f:
        data["RA"] = np.load(f)
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_DEC.npy", "rb") as f:
        data["DEC"] = np.load(f)
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_HPXPIXEL_HPXNSIDE_{NSIDE}.npy", "rb") as f:
        data["HPXPIXEL"] = np.load(f)

    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_BGS_TARGET.npy", "rb") as f:
        bgs_target_bitmasks = np.load(f)

    data["BGS_TARGET"] = np.array([
        BGS_BRIGHT if "BGS_BRIGHT" in bgs_mask.names(i)
        else (BGS_FAINT if "BGS_FAINT" in bgs_mask.names(i)
            else BGS_WISE
        )
        for i in bgs_target_bitmasks
    ])

    fluxes = {}
    mw_transmissions = {}
    mags = {}
    for band in BANDS:
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_FLUX_{band.upper()}.npy", "rb") as f:
            fluxes = np.load(f)
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{region}/targets_MW_TRANSMISSION_{band.upper()}.npy", "rb") as f:
            mw_transmissions = np.load(f)

        data[f"MAG_{band.upper()}"] = np.array(mag_from_flux(fluxes, mw_transmissions, extinction_correction=extinction_correction))

    # Remove spurious objects
    data = remove_spurious_objects(data)

    return data
