from matplotlib import pyplot as plt
import numpy as np

BASS_MzLS = "BASS-MzLS"
DECaLS_NGC = "DECaLS-NGC"
DECaLS_SGC = "DECaLS-SGC"

DESI_REGION = BASS_MzLS

COLORS = {
    14: "C0",
    15: "C1",
    16: "C2",
    17: "C3",
    18: "C4",
    19: "C5",
}

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048

# TODO determine run_id via some better way, or loop through all existing run_id values
run_id = 100

PATH_SHAM_2PCF = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{PARTICLE_COUNT_PINOCCHIO}cubed/{run_id}/2PCF/"
PATH_DESI_LS_2PCF = "/cluster/scratch/lmachado/DataProducts/2PCF/"


# Parameters for plotting
plt.rcParams["font.size"] = "12"
plt.rcParams["figure.figsize"] = (8, 6)
LINEWIDTH = 2

rmag_bins = [
    #[14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

# First plot: compare SHAM 2PCF (all galaxies) vs. DESI LS 2PCF (all galaxies)
# This plot does NOT include any blue vs. red separation
for rmag_low, rmag_high in rmag_bins:
    sham_bins_filename = f"simulated_total_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
    sham_wtheta_filename = f"simulated_total_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

    desi_bins_filename = f"{DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy"
    desi_wtheta_filename = f"{DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy"

    with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
        sham_bins = np.load(f)
    with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
        sham_wtheta = np.load(f)

    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    plt.plot(
        desi_bins,
        desi_wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
        label=f"{rmag_low:.1f} < r < {rmag_high:.1f}",
    )

    plt.plot(
        sham_bins,
        sham_wtheta,
        linewidth=LINEWIDTH,
        linestyle="dashed",
        color=COLORS[int(rmag_low)],
    )

plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$w(\theta)$")
plt.title(f"DESI LS (-) vs. SHAM (--)\n{DESI_REGION}")
plt.xlim([2e-3, 20])
plt.ylim([1e-4, 100])
plt.grid()

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_compareSHAMvsDESI_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()


# --------------------
# Second plot: compare SHAM 2PCF total vs. red vs. blue galaxies
LINESTYLES = {
    "total": "solid",
    "blue": "dashed",
    "red": "dotted",
}

for rmag_low, rmag_high in rmag_bins:
    for color_name, linestyle in LINESTYLES.items():
        sham_bins_filename = f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
            sham_bins = np.load(f)
        with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
            sham_wtheta = np.load(f)

        # Due to 1-halo term issue,
        # only plot BLUE 2PCF
        # above a certain theta
        if color_name == "blue":
            theta_cut = 0.1
            ids = np.where(sham_bins >= theta_cut)[0]
            sham_bins = sham_bins[ids]
            sham_wtheta = sham_wtheta[ids]

        # Only use label for total,
        # so we don't repeat labels
        label = None
        if color_name == "total":
            label = f"{rmag_low:.1f} < r < {rmag_high:.1f}"

        plt.plot(
            sham_bins,
            sham_wtheta,
            linewidth=LINEWIDTH,
            linestyle=linestyle,
            color=COLORS[int(rmag_low)],
            label=label,
        )

plt.legend(loc="upper right")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$w(\theta)$")
plt.title(f"SHAM Comparison: all (-) vs. blue (--) vs. red (...)\n{DESI_REGION}")
plt.xlim([2e-3, 20])
plt.ylim([1e-4, 100])
plt.grid()

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_SHAMRedvsBlue_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()
