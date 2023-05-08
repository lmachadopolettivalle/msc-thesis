from cycler import cycler
from matplotlib import pyplot as plt
import numpy as np

import directories

from manage_parameter_space import get_details_of_run

# Parameters for plotting
plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (8, 15)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"
LINEWIDTH = 2

# Angular scales of trust for each r band range
SCALES_OF_TRUST = {
    (15, 16): (0.2, 1),
    (16, 17): (0.15, 1),
    (17, 18): (0.1, 1),
    (18, 19): (0.09, 1),
    (19, 19.5): (0.08, 1),
}

color_cycler = cycler(color=[
    "#1965B0",
    "#F1932D",
    "#4EB265",
    "#DC050C",
    "#AA6F9E",
])

# Choose whether to fix mlimit or tquench. Plots will fix this variable, and vary the other one
# run_ids have been found manually to show the most variation with this fixed parameter
MLIMIT = "mass_cut"
TQUENCH = "quenching_time"

FIXED_PARAMETER = TQUENCH # Change between MLIMIT and TQUENCH to obtain different figures

if FIXED_PARAMETER == MLIMIT:
    FIXED_VALUE = 8e12
    RUN_IDS = [142, 155, 143, 151]
    VARIED_PARAMETER = TQUENCH
    VARIED_PARAMETER_LABEL = "$t_{{quench}}$"

elif FIXED_PARAMETER == TQUENCH:
    FIXED_VALUE = 1.5
    RUN_IDS = [152, 154, 147, 155]
    VARIED_PARAMETER = MLIMIT
    VARIED_PARAMETER_LABEL = "$M_{{limit}}$"
else:
    raise ValueError("Invalid setting for the fixed parameter.")


DESI_REGION = directories.FULLSKY

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

rmag_bins = [
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

PATH_DESI_LS_2PCF = "/cluster/scratch/lmachado/DataProducts/2PCF/"

# Create blank figure for plots
fig, axs = plt.subplots(nrows=len(rmag_bins), ncols=1, sharex=True)
#fig.subplots_adjust(wspace=0.1, hspace=0.15)

for ax, (rmag_low, rmag_high) in zip(axs, rmag_bins):
    ax.set_prop_cycle(color_cycler)

    # Obtain data from DESI LS
    desi_bins_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}.npy"
    desi_wtheta_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}.npy"

    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    ax.plot(
        desi_bins,
        [0] * len(desi_bins),
        linewidth=LINEWIDTH,
        color="black",
    )

    # Plot ratio between SHAM and DESI LS 2PCF for each run_id
    for run_id in RUN_IDS:

        # Obtain SHAM data and plot its ratio to DESI LS 2PCF
        run_details = get_details_of_run(run_id)

        num_mass_bins = run_details["num_mass_bins"]
        num_z_bins = run_details["num_z_bins"]
        mlimit = run_details["mass_cut"]
        tquench = run_details["quenching_time"]

        PATH_SHAM_2PCF = directories.path_2PCF(
            particle_count=PARTICLE_COUNT_PINOCCHIO,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=DESI_REGION,
            run_id=run_id,
        )

        sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
            sham_bins = np.load(f)
        with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
            sham_wtheta = np.load(f)

        trusted_indices = np.where(
            (sham_bins >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) & (sham_bins <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
        )[0]

        ax.plot(
            desi_bins[trusted_indices],
            sham_wtheta[trusted_indices] / desi_wtheta[trusted_indices] - 1,
            linewidth=LINEWIDTH,
            linestyle="solid",
            label=(VARIED_PARAMETER_LABEL + f" = {run_details[VARIED_PARAMETER]:.1e} $M_{{\odot}}/h$") if VARIED_PARAMETER == MLIMIT else (VARIED_PARAMETER_LABEL + f" = {run_details[VARIED_PARAMETER]:.1f} Gyr"),
        )

        ax.text(
            0.6, # x coord
            -0.37, # y coord
            f"{rmag_low:.1f} < r < {rmag_high:.1f}",
            bbox=dict(facecolor="none", edgecolor="black", pad=6),
            ha="center",
            va="center",
        )

    ax.set_xscale("log")

    ax.set_xlim([8e-2, 1])
    ax.set_ylim([-0.5, 0.5])

    ax.grid()

axs[0].legend(
    loc="upper left",
    bbox_to_anchor=(1.04, 1.0),
)

# Label for the fixed variable
if FIXED_PARAMETER == MLIMIT:
    fixed_value_label = f"$M_{{limit}}$ = {mlimit:.1e} $M_{{\odot}}/h$" 
else:
    fixed_value_label = f"$t_{{quench}}$ = {tquench:.1f} Gyr" 


axs[0].text(
    1.2, # x coord
    -0.7, # y coord
    f"$m_{{bins}}$ = {num_mass_bins}\n$z_{{bins}}$ = {num_z_bins}\n" + fixed_value_label,
    bbox=dict(facecolor="none", edgecolor="black", pad=6),
    ha="left",
    va="center",
)

axs[-1].set_xlabel(r"$\theta$ [deg]")
axs[2].set_ylabel(r"$w(\theta)_{SHAM}$ / $w(\theta)_{data} - 1$", labelpad=8)

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_impactSHAM_fixed{FIXED_PARAMETER}_{DESI_REGION}.pdf")

plt.show()
