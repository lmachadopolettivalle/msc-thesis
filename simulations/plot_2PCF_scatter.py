from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate

import directories

from manage_parameter_space import get_details_of_run

DESI_REGION = directories.FULLSKY
PINOCCHIO_REGION = "fullsky"

# Whether to show narrow or wide theta range in 2PCF comparison
NARROW_THETA_RANGE = False

# Whether to use r-primed
# Only applies to BASS_MzLS objects
USE_MAG_R = True


if DESI_REGION == directories.FULLSKY:
    EQUIVALENT_DESI_REGION = (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC)
else:
    EQUIVALENT_DESI_REGION = DESI_REGION

REGION_LINESTYLE = {
    directories.BASS_MzLS: "solid",
    directories.DECaLS_NGC: "dashed",
    directories.DECaLS_SGC: "dashdot",
    directories.FULLSKY: "solid",
}
REGION_MARKERS = {
    directories.BASS_MzLS: None,
    directories.DECaLS_NGC: "o",
    directories.DECaLS_SGC: "^",
    directories.FULLSKY: "*",
}

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5

# Select run_id for 2PCF visualization
run_id = 148
run_details = get_details_of_run(run_id)

PATH_DESI_LS_2PCF = "/cluster/scratch/lmachado/DataProducts/2PCF/"

# Parameters for plotting
plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (8, 9)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"
LINEWIDTH = 2

rmag_bins = [
    #[14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

COLORS = {
    15: "#1965B0", # blue brighter
    16: "#F1932D", # orange
    17: "#4EB265", # green brighter
    18: "#DC050C", # red brighter
    19: "#AA6F9E", # purple a bit brighter
}

# Angular scales of trust for each r band range
SCALES_OF_TRUST = {
    (15, 16): (0.22, 3),
    (16, 17): (0.15, 3),
    (17, 18): (0.1, 3),
    (18, 19): (0.09, 3),
    (19, 19.5): (0.08, 3),
}

fig, ax = plt.subplots(1, 1)

seeds = [None] + [f"{i * 111111}" for i in range(1, 6+1)]

for rmag_low, rmag_high in rmag_bins:
    desi_bins_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}_noextinction.npy"
    desi_wtheta_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}_noextinction.npy"
    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    ax.plot(
        desi_bins,
        desi_wtheta,
        label=f"DESI LS, {rmag_low:.1f} < r < {rmag_high:.1f}",
        lw=2,
    )

    sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
    sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

    sham_wtheta_values = []

    for seed in seeds:
        PATH_SHAM_2PCF = directories.path_2PCF(
            particle_count=PARTICLE_COUNT_PINOCCHIO,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=DESI_REGION,
            run_id=run_id,
            seed=seed,
        )

        with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
            sham_bins = np.load(f)
        with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
            sham_wtheta = np.load(f)

        sham_wtheta_values.append(sham_wtheta)

    sham_wtheta_values = np.array(sham_wtheta_values)
    sham_wtheta_mean = sham_wtheta_values.mean(axis=0)
    sham_wtheta_std = sham_wtheta_values.std(axis=0)

    ax.fill_between(
        sham_bins,
        sham_wtheta_mean - sham_wtheta_std,
        sham_wtheta_mean + sham_wtheta_std,
        alpha=0.2,
        lw=2,
        label=f"PINOCCHIO, {rmag_low:.1f} < r < {rmag_high:.1f}"
    )

ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\theta$ [deg]")
ax.set_ylabel(r"$w(\theta)$")
ax.set_xlim([0.08, 3])
ax.set_ylim([1e-5, 2])
plt.show()
