# Read file with data from selected targets
import numpy as np
from desitarget.targets import bgs_mask
from mask import mask
from constants import *

# Parameters for magnitude computation
CLIP_FLUX = 1e-16

def mag_from_flux(fluxes, mw_transmissions, extinction_correction=True):
    if extinction_correction:
        mags = 22.5 - 2.5 * np.log10((fluxes / mw_transmissions).clip(CLIP_FLUX))
    else:
        mags = 22.5 - 2.5 * np.log10(fluxes.clip(CLIP_FLUX))

    return mags

# Equation 2 in Zarrouk+21
def rmag_bass_correction(unprimed_rmags, unprimed_gmags):
    assert len(unprimed_rmags) == len(unprimed_gmags)
    return unprimed_rmags - 0.039 * (unprimed_gmags - unprimed_rmags) + 0.011

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

def load_processed_target_data(regions=ALL_REGIONS, extinction_correction=True, apply_mask=True):
    # Return dict, with values being arrays of the following fields for targets:
    # RA, DEC, Magnitudes, BGS Bright or Faint
    NSIDE = 512

    data = {
        "RA": np.array([]),
        "DEC": np.array([]),
        "HPXPIXEL": np.array([], dtype=int),
        "BGS_TARGET": np.array([]),
        "MAG_R": np.array([]),
        "MAG_R_PRIMED": np.array([]), # For BASS r-mag, apply correction via equation 2 from Zarrouk+21
        "MAG_G": np.array([]),
        "MAG_Z": np.array([]),
        "MORPHTYPE": np.array([]),
    }

    sky_regions = set()
    for region in regions:
        sky_regions.add(map_region_to_north_south(region))

    for sky_region in sky_regions:
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_RA.npy", "rb") as f:
            data["RA"] = np.concatenate((data["RA"], np.load(f)))
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_DEC.npy", "rb") as f:
            data["DEC"] = np.concatenate((data["DEC"], np.load(f)))
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_HPXPIXEL_HPXNSIDE_{NSIDE}.npy", "rb") as f:
            data["HPXPIXEL"] = np.concatenate((data["HPXPIXEL"], np.load(f)), dtype=int)
        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_MORPHTYPE.npy", "rb") as f:
            data["MORPHTYPE"] = np.concatenate((data["MORPHTYPE"], np.load(f)))

        with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_BGS_TARGET.npy", "rb") as f:
            bgs_target_bitmasks = np.load(f)

        bgs_targets_array = np.array([
            BGS_BRIGHT if "BGS_BRIGHT" in bgs_mask.names(i)
            else (BGS_FAINT if "BGS_FAINT" in bgs_mask.names(i)
                else BGS_WISE
            )
            for i in bgs_target_bitmasks
        ])

        data["BGS_TARGET"] = np.concatenate((data["BGS_TARGET"], bgs_targets_array))

        fluxes = {}
        mw_transmissions = {}
        mags = {}
        for band in sorted(BANDS):
            # BANDS must be sorted, to ensure we have computed g-band info before r-band,
            # since g-band is required for r-mag correction from Zarrouk+21

            with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_FLUX_{band.upper()}.npy", "rb") as f:
                fluxes = np.load(f)
            with open(f"/cluster/scratch/lmachado/DataProducts/targets/{sky_region}/targets_MW_TRANSMISSION_{band.upper()}.npy", "rb") as f:
                mw_transmissions = np.load(f)

            mag_array = np.array(mag_from_flux(fluxes, mw_transmissions, extinction_correction=extinction_correction))
            data[f"MAG_{band.upper()}"] = np.concatenate((data[f"MAG_{band.upper()}"], mag_array))

            # Apply correction to r-mag from Zarrouk+21
            # Only apply correction to northern targets,
            # since correction is meant to make BASS r-mag values match
            # more closely with DECaLS r-mag values
            if band.lower() == "r":
                if sky_region == "north":
                    # To get corresponding g-mag values, read from data["MAG_G"],
                    # knowing that the latest values have been appended for this same sky region,
                    # with the same number of targets
                    primed_rmag_array = rmag_bass_correction(
                        mag_array,
                        data["MAG_G"][-1*len(mag_array):]
                    )
                    data["MAG_R_PRIMED"] = np.concatenate((data["MAG_R_PRIMED"], primed_rmag_array))
                elif sky_region == "south":
                    data["MAG_R_PRIMED"] = np.concatenate((data["MAG_R_PRIMED"], mag_array))
                else:
                    raise ValueError(f"Sky region cannot be {sky_region}, only north or south!")

    # Remove spurious objects
    data = remove_spurious_objects(data)

    # Convert MORPHTYPE from bytes to str
    data["MORPHTYPE"] = np.vectorize(lambda x: x.decode("utf-8"))(data["MORPHTYPE"])

    # Remove objects of certain MORPHTYPE (e.g. PSF, which are stars)
    # NOTE we invert the mask, to keep only the other MORPHTYPE values
    good_morphtypes_mask = np.isin(data["MORPHTYPE"], ["PSF"], invert=True)
    data = {
        k: v[good_morphtypes_mask]
        for k, v in data.items()
    }

    # If requested, apply mask to targets
    if apply_mask:
        targets_masked = mask(data["HPXPIXEL"], NSIDE, regions=regions)
        targets_ids_in_mask = np.where(targets_masked > 0)[0]

        data = {
            k: v[targets_ids_in_mask]
            for k, v in data.items()
        }

    return data
