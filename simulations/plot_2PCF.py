from matplotlib import pyplot as plt
import numpy as np

import directories

from manage_parameter_space import get_details_of_run

# Whether to use r-primed
# Only applies to BASS_MzLS objects
USE_MAG_R = True

DESI_REGION = directories.DECaLS_SGC

# Automaticaly set USE_MAG_R to False,
# if not in BASS-MzLS region
if DESI_REGION != directories.BASS_MzLS:
    USE_MAG_R = False

REGION_LINESTYLE = {
    directories.BASS_MzLS: "solid",
    directories.DECaLS_NGC: "dashed",
    directories.DECaLS_SGC: "dashdot",
    (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC): "solid",
}
REGION_MARKERS = {
    directories.BASS_MzLS: None,
    directories.DECaLS_NGC: "o",
    directories.DECaLS_SGC: "^",
    (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC): "*",
}

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

# Select run_id for 2PCF visualization
run_id = 146
run_details = get_details_of_run(run_id)
# Optionally, set a second run_id for a comparison between the two runs
# If not wanted, set the baseline_run_id to None
baseline_run_id = None

PATH_SHAM_2PCF = directories.path_2PCF(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_REGION,
    run_id=run_id,
)
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

COLORS = {
    #14: "C0",
    15: "C0",
    16: "C1",
    17: "C2",
    18: "C3",
    19: "C4",
}


# First plot: compare SHAM 2PCF (all galaxies) vs. DESI LS 2PCF (all galaxies)
# This plot does NOT include any blue vs. red separation
fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

for rmag_low, rmag_high in rmag_bins:
    sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
    sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

    desi_bins_filename = f"{DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R else 'unprimed'}.npy"
    desi_wtheta_filename = f"{DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R else 'unprimed'}.npy"

    with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
        sham_bins = np.load(f)
    with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
        sham_wtheta = np.load(f)

    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    ax_top.plot(
        desi_bins,
        desi_wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
        label=f"{rmag_low:.1f} < r < {rmag_high:.1f}",
    )

    ax_top.plot(
        sham_bins,
        sham_wtheta,
        linewidth=LINEWIDTH,
        linestyle="dashed",
        color=COLORS[int(rmag_low)],
    )

    ax_bottom.plot(
        desi_bins,
        sham_wtheta / desi_wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
    )
    ax_bottom.plot(
        desi_bins,
        [1] * len(desi_bins),
        linewidth=LINEWIDTH,
        color="black",
    )

ax_top.legend()

ax_top.set_xscale("log")
ax_top.set_yscale("log")
ax_bottom.set_xscale("log")

ax_top.set_ylabel(r"$w(\theta)$")
ax_bottom.set_xlabel(r"$\theta$ [deg]")
ax_bottom.set_ylabel("Ratio")

ax_top.set_xlim([1e-2, 5])
ax_top.set_ylim([1e-4, 10])
ax_bottom.set_ylim([0.5, 1.5])

ax_top.grid()
ax_bottom.grid()

plt.subplots_adjust(hspace=0.1)

fig.suptitle(f"DESI LS (-) vs. SHAM (--)\n{DESI_REGION}")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_compareSHAMvsDESI_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

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
        sham_bins_filename = f"simulated_{color_name}_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_{color_name}_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        try:
            with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
                sham_bins = np.load(f)
            with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
                sham_wtheta = np.load(f)
        except FileNotFoundError as e:
            print(e)
            print(f"Skipping color_name {color_name}...")
            break

        # Due to 1-halo term issue,
        # only plot BLUE 2PCF
        # above a certain theta
        """
        if color_name == "blue":
            theta_cut = 0.1
            ids = np.where(sham_bins >= theta_cut)[0]
            sham_bins = sham_bins[ids]
            sham_wtheta = sham_wtheta[ids]
        """

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
plt.xlim([1e-2, 5])
plt.ylim([1e-4, 10])
plt.grid()

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_SHAMRedvsBlue_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()


# --------------------
# Third plot: compare SHAM for each of the three regions (BASS-MzLS, DECaLS-NGC, DECaLS-SGC)
# This plot does NOT include any blue vs. red separation
fig, ax = plt.subplots(1, 1, figsize=(14, 8))

for region in (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC):
    # Reset colors, to use same colors for any given rmag range across all regions
    ax.set_prop_cycle(None)

    for rmag_low, rmag_high in rmag_bins:
        path_2PCF_region = directories.path_2PCF(
            particle_count=PARTICLE_COUNT_PINOCCHIO,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=region,
            run_id=run_id,
        )
        sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if (USE_MAG_R and (region==directories.BASS_MzLS)) else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if (USE_MAG_R and (region==directories.BASS_MzLS)) else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{path_2PCF_region}/{sham_bins_filename}", "rb") as f:
            bins = np.load(f)
        with open(f"{path_2PCF_region}/{sham_wtheta_filename}", "rb") as f:
            wtheta = np.load(f)

        ax.plot(
            bins,
            wtheta,
            linewidth=LINEWIDTH,
            linestyle=REGION_LINESTYLE[region],
            color=COLORS[int(rmag_low)],
            marker=REGION_MARKERS[region],
            label=f"{rmag_low:.1f} < r < {rmag_high:.1f}, {region}",
        )

ax.legend(loc="lower left", ncol=3)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$\theta$ [deg]")
ax.set_ylabel(r"$w(\theta)$")

ax.set_xlim([1e-2, 5])
ax.set_ylim([1e-4, 10])

ax.grid()

ax.set_title("SHAM Comparison between Regions")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_compare_SHAMregions_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()

# --------------------
# Fourth plot: Compare r-primed and r 2PCF for BASS
# This plot does NOT include any blue vs. red separation
fig, ax = plt.subplots(1, 1)

for is_prime in (True, False):
    for rmag_low, rmag_high in rmag_bins:
        path_2PCF_region = directories.path_2PCF(
            particle_count=PARTICLE_COUNT_PINOCCHIO,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=directories.BASS_MzLS,
            run_id=run_id,
        )
        sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if is_prime else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if is_prime else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{path_2PCF_region}/{sham_bins_filename}", "rb") as f:
            bins = np.load(f)
        with open(f"{path_2PCF_region}/{sham_wtheta_filename}", "rb") as f:
            wtheta = np.load(f)

        ax.plot(
            bins,
            wtheta,
            linewidth=LINEWIDTH,
            linestyle=("solid" if is_prime else "dashed"),
            color=COLORS[int(rmag_low)],
            label=f"{rmag_low:.1f} < r < {rmag_high:.1f}, {'Primed' if is_prime else 'Unprimed'}",
        )

ax.legend(loc="lower left", ncol=2)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$\theta$ [deg]")
ax.set_ylabel(r"$w(\theta)$")

ax.set_xlim([1e-2, 5])
ax.set_ylim([1e-4, 10])

ax.grid()

ax.set_title("SHAM 2PCF, r-primed vs. r, BASS-MzLS")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_compare_rprimed_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()

# --------------------
# Fifth plot: Compare 2PCF between two run IDs
if baseline_run_id is None:
    print("Baseline run ID has not been set. Exiting program.")
    exit()

baseline_run_details = get_details_of_run(baseline_run_id)

PATH_BASELINE_SHAM_2PCF = directories.path_2PCF(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_REGION,
    run_id=baseline_run_id,
)

fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

for rmag_low, rmag_high in rmag_bins:
    sham_bins_filename = f"simulated_total_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
    sham_wtheta_filename = f"simulated_total_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

    with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
        bins = np.load(f)
    with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
        wtheta = np.load(f)

    with open(f"{PATH_BASELINE_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
        baseline_bins = np.load(f)
    with open(f"{PATH_BASELINE_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
        baseline_wtheta = np.load(f)

    ax_top.plot(
        bins,
        wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
        label=f"{rmag_low:.1f} < r < {rmag_high:.1f}",
    )

    ax_top.plot(
        baseline_bins,
        baseline_wtheta,
        linewidth=LINEWIDTH,
        linestyle="dashed",
        color=COLORS[int(rmag_low)],
    )

    ax_bottom.plot(
        bins,
        wtheta / baseline_wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
    )
    ax_bottom.plot(
        bins,
        [1] * len(bins),
        linewidth=LINEWIDTH,
        color="black",
    )

ax_top.legend()

ax_top.set_xscale("log")
ax_top.set_yscale("log")
ax_bottom.set_xscale("log")

ax_top.set_ylabel(r"$w(\theta)$")
ax_bottom.set_xlabel(r"$\theta$ [deg]")
ax_bottom.set_ylabel("Ratio")

ax_top.set_xlim([1e-2, 5])
ax_top.set_ylim([1e-4, 10])

ax_top.grid()
ax_bottom.grid()

plt.subplots_adjust(hspace=0.1)

fig.suptitle(f"Solid: m_bins = {run_details['num_mass_bins']}, z_bins = {run_details['num_z_bins']}\nDashed: m_bins = {baseline_run_details['num_mass_bins']}, z_bins = {baseline_run_details['num_z_bins']}")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_compareRUNS_{PARTICLE_COUNT_PINOCCHIO}_{run_id}_{baseline_run_id}.pdf")

plt.show()
