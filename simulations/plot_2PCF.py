from matplotlib import pyplot as plt
import numpy as np

import directories

from manage_parameter_space import get_details_of_run

DESI_REGION = directories.FULLSKY

# Whether to show narrow or wide theta range in 2PCF comparison
NARROW_THETA_RANGE = True

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
PINOCCHIO_REGION = "fullsky"

# Select run_id for 2PCF visualization
run_id = 154
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
    #14: "C0",
    15: "#E8601C",
    16: "#F6C141",
    17: "#90C987",
    18: "#5289C7",
    19: "#AE76A3",
}

# Angular scales of trust for each r band range
SCALES_OF_TRUST = {
    (15, 16): (0.22, 3),
    (16, 17): (0.15, 3),
    (17, 18): (0.1, 3),
    (18, 19): (0.09, 3),
    (19, 19.5): (0.08, 3),
}

# --------------------
# First plot: compare SHAM 2PCF (all galaxies) vs. DESI LS 2PCF (all galaxies)
# This plot does NOT include any blue vs. red separation
fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})

plotted_lines = []

for rmag_low, rmag_high in rmag_bins:
    sham_bins_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
    sham_wtheta_filename = f"simulated_total_2PCF_{'rprimed_' if USE_MAG_R else ''}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

    desi_bins_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}.npy"
    desi_wtheta_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if USE_MAG_R and (DESI_REGION in (directories.BASS_MzLS, directories.FULLSKY)) else 'unprimed'}.npy"

    with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
        sham_bins = np.load(f)
    with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
        sham_wtheta = np.load(f)

    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    # If focusing on narrow theta range, filter SHAM plots to the scales of trust
    if NARROW_THETA_RANGE:
        trusted_indices = np.where((sham_bins >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) & (sham_bins <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1]))[0]
    else:
        trusted_indices = list(range(len(sham_bins)))

    desi_line, = ax_top.plot(
        desi_bins[trusted_indices],
        desi_wtheta[trusted_indices],
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
    )

    sham_line, = ax_top.plot(
        sham_bins[trusted_indices],
        sham_wtheta[trusted_indices],
        linewidth=LINEWIDTH,
        linestyle="dashed",
        color=COLORS[int(rmag_low)],
    )

    plotted_lines.append([desi_line, sham_line])

    ax_bottom.plot(
        desi_bins[trusted_indices],
        sham_wtheta[trusted_indices] / desi_wtheta[trusted_indices] - 1,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
    )
    ax_bottom.plot(
        desi_bins,
        [0] * len(desi_bins),
        linewidth=LINEWIDTH,
        color="black",
    )

ax_top.set_xscale("log")
ax_top.set_yscale("log")
ax_bottom.set_xscale("log")

ax_top.set_ylabel(r"$w(\theta)$")
ax_bottom.set_xlabel(r"$\theta$ [deg]")
ax_bottom.set_ylabel("Fractional Difference")


if NARROW_THETA_RANGE:
    ax_top.set_xlim([1e-1, 1])
    ax_top.set_ylim([1e-3, 10])
    ax_bottom.set_ylim([-0.2, 0.2])
    rmags_legend_loc = "upper left"
    ax_bottom.set_xticks(
        [0.1, 0.2, 0.3, 0.4, 0.6, 1],
        ["0.1", "0.2", "0.3", "0.4", "0.6", "1.0"],
    )
else:
    ax_top.set_xlim([1e-2, 20])
    ax_top.set_ylim([1e-4, 10])
    ax_bottom.set_ylim([-1, 1])
    rmags_legend_loc = "upper right"

rmags_legend = ax_top.legend(
    [i[0] for i in plotted_lines],
    [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in rmag_bins],
    loc=rmags_legend_loc,
)

_ = ax_top.legend(
    plotted_lines[0],
    ["DESI LS", "SHAM"],
    loc="lower left",
)

ax_top.add_artist(rmags_legend)


ax_top.grid()
ax_bottom.grid()

plt.subplots_adjust(hspace=0.1)

#fig.suptitle(f"DESI LS (-) vs. SHAM (--)\n{DESI_REGION}")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_compareSHAMvsDESI_{PARTICLE_COUNT_PINOCCHIO}_{run_id}{'narrow' if NARROW_THETA_RANGE else ''}.pdf")

plt.show()

exit() # TODO

# --------------------
# Second plot: compare SHAM for each of the four regions (BASS-MzLS, DECaLS-NGC, DECaLS-SGC)
# This plot does NOT include any blue vs. red separation
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

plotted_lines = []

for region in (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC, directories.FULLSKY):
    # Reset colors, to use same colors for any given rmag range across all regions
    ax.set_prop_cycle(None)

    tmp_plotted_lines = []

    for rmag_low, rmag_high in rmag_bins:
        path_2PCF_region = directories.path_2PCF(
            particle_count=PARTICLE_COUNT_PINOCCHIO,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=region,
            run_id=run_id,
        )
        if USE_MAG_R:
            primed_text = "rprimed_"
        else:
            primed_text = ""
        sham_bins_filename = f"simulated_total_2PCF_{primed_text}{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_{primed_text}{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{path_2PCF_region}/{sham_bins_filename}", "rb") as f:
            bins = np.load(f)
        with open(f"{path_2PCF_region}/{sham_wtheta_filename}", "rb") as f:
            wtheta = np.load(f)

        line, = ax.plot(
            bins,
            wtheta,
            linewidth=LINEWIDTH,
            linestyle=REGION_LINESTYLE[region],
            color=COLORS[int(rmag_low)],
            marker=REGION_MARKERS[region],
            label=f"{rmag_low:.1f} < r < {rmag_high:.1f}, {region}",
        )

        tmp_plotted_lines.append(line)

    plotted_lines.append(tmp_plotted_lines)

rmags_legend = ax.legend(
    plotted_lines[0],
    [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in rmag_bins],
    loc="upper right",
)

_ = ax.legend(
    [i[0] for i in plotted_lines],
    [directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC, directories.FULLSKY],
    loc="lower left",
)

ax.add_artist(rmags_legend)

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
# Third plot: Compare r-primed and r 2PCF for BASS
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

ax.set_xlim([1e-2, 20])
ax.set_ylim([1e-4, 10])

ax.grid()

ax.set_title("SHAM 2PCF, r-primed vs. r, BASS-MzLS")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_compare_rprimed_{PARTICLE_COUNT_PINOCCHIO}_{run_id}.pdf")

plt.show()

# --------------------
# Fourth plot: Compare 2PCF between two run IDs
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

plotted_lines = []

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

    first_line, = ax_top.plot(
        bins,
        wtheta,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
        #label=f"{rmag_low:.1f} < r < {rmag_high:.1f}",
    )

    baseline_line, = ax_top.plot(
        baseline_bins,
        baseline_wtheta,
        linewidth=LINEWIDTH,
        linestyle="dashed",
        color=COLORS[int(rmag_low)],
    )

    plotted_lines.append([first_line, baseline_line])

    ax_bottom.plot(
        bins,
        wtheta / baseline_wtheta - 1,
        linewidth=LINEWIDTH,
        linestyle="solid",
        color=COLORS[int(rmag_low)],
    )
    ax_bottom.plot(
        bins,
        [0] * len(bins),
        linewidth=LINEWIDTH,
        color="black",
    )

rmags_legend = ax_top.legend(
    [i[0] for i in plotted_lines],
    [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in rmag_bins],
    loc="upper right",
)

_ = ax_top.legend(
    plotted_lines[0],
    [
        f"$m_{{bins}}$ = {run_details['num_mass_bins']}, $z_{{bins}}$ = {run_details['num_z_bins']}",
        f"$m_{{bins}}$ = {baseline_run_details['num_mass_bins']}, $z_{{bins}}$ = {baseline_run_details['num_z_bins']}",
    ],
    loc="lower left",
)

ax_top.add_artist(rmags_legend)

ax_top.set_xscale("log")
ax_top.set_yscale("log")
ax_bottom.set_xscale("log")

ax_top.set_ylabel(r"$w(\theta)$")
ax_bottom.set_xlabel(r"$\theta$ [deg]")
ax_bottom.set_ylabel("Fractional Difference")

ax_top.set_xlim([0.1, 1])
ax_top.set_ylim([1e-3, 10])
ax_bottom.set_ylim([-0.2, 0.2])

ax_bottom.set_xticks(
    [0.1, 0.2, 0.3, 0.4, 0.6, 1],
    ["0.1", "0.2", "0.3", "0.4", "0.6", "1.0"],
)

ax_top.grid()
ax_bottom.grid()

plt.subplots_adjust(hspace=0.1)

#fig.suptitle(f"Solid: m_bins = {run_details['num_mass_bins']}, z_bins = {run_details['num_z_bins']}\nDashed: m_bins = {baseline_run_details['num_mass_bins']}, z_bins = {baseline_run_details['num_z_bins']}")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_{DESI_REGION}_compareRUNS_{PARTICLE_COUNT_PINOCCHIO}_{run_id}_{baseline_run_id}.pdf")

plt.show()
