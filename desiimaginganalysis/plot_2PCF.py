from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

from constants import *

# Plots generated:
# r-primed vs unprimed, for BASS only, and for BASS + DECaLS-NGC + DECaLS-SGC
# overplotted curves for the 3 regions, without gamma factor
# overplotted curves for the 3 regions, with gamma factor, and overplot Zarrouk
# overplot Zarrouk and BASS + DECaLS-NGC + DECaLS-SGC, using gammas

# Decide if show every plot or just save and clear figure
SHOW = True

# Parameters for plotting
plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"
LINEWIDTH = 2

BRIGHT = "Bright"
FAINT = "Faint"

LINESTYLE = {
    BASS_MzLS: "solid",
    DECaLS_NGC: "dashed",
    DECaLS_SGC: "dashdot",
    ALL_REGIONS: "solid",
}

MARKERS = {
    BASS_MzLS: None,
    DECaLS_NGC: "o",
    DECaLS_SGC: "^",
    ALL_REGIONS: "*",
}

COLORS = {
    15: "#1965B0", # blue brighter
    16: "#F1932D", # orange
    17: "#4EB265", # green brighter
    18: "#DC050C", # red brighter
    19: "#AA6F9E", # purple a bit brighter
}

RMAG_BINS = [
    (15, 16),
    (16, 17),
    (17, 18),
    (18, 19),
    (19, 19.5),
]

# Obtain list of 2PCF files available
PATH_TO_2PCF_FILES = "/cluster/scratch/lmachado/DataProducts/2PCF/"

# Focus on desired region, on _bins_ files, and on bright vs. faint
FILENAMES = {
    region: [
        f for f in sorted(os.listdir(PATH_TO_2PCF_FILES), key=lambda x:x[::-1])
        if f.startswith(str(region)) and ("_bins_" in f)
    ]
    
    for region in [
        BASS_MzLS,
        DECaLS_NGC,
        DECaLS_SGC,
        ALL_REGIONS,
    ]
}

TYPE_FILENAMES = {
    brightness: {
        region: [
            f for f in filenames
            if brightness in f
        ]

        for region, filenames in FILENAMES.items()
    }
    for brightness in (BRIGHT, FAINT)
}

def fit_power_law(xs, ys):
    def func_powerlaw(x, exp, norm):
        return x**exp * norm

    # Only fit using 1e-3 < theta < 1
    ids = np.where((1e-3 < ys) & (ys < 1))[0]

    popt, pcov = curve_fit(
        func_powerlaw,
        xs[ids],
        ys[ids],
        p0=np.asarray([1.6, 0.05]), # Initial guess for parameters
        maxfev=2000,
    )

    return popt, pcov

def get_rmags_isprimed_from_filename(filename):
    if "unprimed" in filename:
        primed_filename_extension = "_unprimed"
    else:
        primed_filename_extension = "_primed"

    rmag_range = filename[
        filename.find("range")+len("range"):
        filename.find(f"{primed_filename_extension}.npy")
    ]
    rmag_low, rmag_high = rmag_range.split("-")

    return (
        float(rmag_low),
        float(rmag_high),
        (primed_filename_extension == "_primed"),
    )


# Gammas from Table 2 of Zarrouk+2021
GAMMAS = {
    BASS_MzLS: {
        15: 1.642,
        16: 1.744,
        17: 1.776,
        18: 1.750,
        19: 1.725,
    },
    DECaLS_NGC: {
        15: 1.761,
        16: 1.761,
        17: 1.746,
        18: 1.742,
        19: 1.720,
    },
    DECaLS_SGC: {
        15: 1.698,
        16: 1.715,
        17: 1.793,
        18: 1.745,
        19: 1.740,
    },
    ALL_REGIONS: {
        # NOTE: the paper does not have gammas per r-range for the entire target set,
        # so I am using the DECaLS_SGC gammas, since it is the region with the most targets
        # Could also consider using an average of the gammas for each of the 3 regions instead.
        15: 1.698,
        16: 1.735,
        17: 1.793,
        18: 1.745,
        19: 1.740,
    },
}



"""
# r-primed vs unprimed for BASS and for all-sky
print("r-primed vs unprimed")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    for region in [
        BASS_MzLS,
        ALL_REGIONS,
    ]:
        filelist = filedicts[region]

        fig, ax = plt.subplots(1, 1)

        for filename in filelist:
            with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
                bins = np.load(f)

            rmag_low, rmag_high, is_primed = get_rmags_isprimed_from_filename(filename)
            if is_primed:
                primed_filename_extension = "_primed"
            else:
                primed_filename_extension = "_unprimed"

            with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}{primed_filename_extension}.npy", "rb") as f:
                wtheta = np.load(f)

            ax.plot(
                bins,
                wtheta,
                linewidth=LINEWIDTH,
                linestyle=("solid" if is_primed else "dashed"),
                color=COLORS[int(rmag_low)],
                label=f"{rmag_low} < r < {rmag_high}, {'primed' if is_primed else 'unprimed'}"
            )

        ax.set_xlim([1e-2, 20])
        ax.set_ylim([1e-4, 10])
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_ylabel(r"$w(\theta)$")
        ax.set_title(f"BGS {brightness} Targets, {region}\nComparison r-primed vs. r-unprimed")

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.legend(loc="upper right", ncol=2)

        ax.grid()

        fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{region}_2PCF_{brightness}_rprimed_vs_unprimed.pdf")

        if SHOW:
            plt.show()
        else:
            plt.clf()
"""


# Overplotted curves for the 3 regions, without gamma factor
print("Compare 3 regions without gamma")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    plotted_lines = []

    fig, ax = plt.subplots(1, 1)

    for region in [
        BASS_MzLS,
        DECaLS_NGC,
        DECaLS_SGC,
    ]:
        filelist = filedicts[region]

        tmp_plotted_lines = []

        # Reset colors, to use same colors for any given rmag range across all regions
        ax.set_prop_cycle(None)

        for filename in filelist:
            if "noextinction" in filename:
                continue
            with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
                bins = np.load(f)

            rmag_low, rmag_high, is_primed = get_rmags_isprimed_from_filename(filename)
            if is_primed:
                continue

            with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
                wtheta = np.load(f)

            line, = ax.plot(
                bins,
                wtheta,
                linewidth=LINEWIDTH,
                linestyle=LINESTYLE[region],
                color=COLORS[rmag_low],
                marker=MARKERS[region],
                #label=f"{rmag_low} < r < {rmag_high}",
            )
            tmp_plotted_lines.append(line)
        plotted_lines.append(tmp_plotted_lines)

    ax.set_xlim([1e-2, 20])
    ax.set_ylim([1e-4, 10])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$")
    #ax.set_title(f"BGS {brightness} Targets\nComparison between Regions (unprimed)")

    ax.set_xscale("log")
    ax.set_yscale("log")

    rmags_legend = ax.legend(
        plotted_lines[0],
        [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in RMAG_BINS],
        loc="upper right",
    )

    _ = ax.legend(
        [i[0] for i in plotted_lines],
        [BASS_MzLS, DECaLS_NGC, DECaLS_SGC],
        loc="lower left",
    )

    ax.add_artist(rmags_legend)

    ax.grid()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{brightness}_region_comparison.pdf")

    if SHOW:
        plt.show()
    else:
        plt.clf()

# overplotted curves for the 3 regions, with gamma factor, and overplot Zarrouk
print("Compare 3 regions with gamma")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    fig, ax = plt.subplots(1, 1)

    for region in [
        BASS_MzLS,
        DECaLS_NGC,
        DECaLS_SGC,
    ]:
        filelist = filedicts[region]

        # Reset colors, to use same colors for any given rmag range across all regions
        ax.set_prop_cycle(None)

        for filename in filelist:
            if "noextinction" in filename:
                continue
            with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
                bins = np.load(f)

            rmag_low, rmag_high, is_primed = get_rmags_isprimed_from_filename(filename)
            if is_primed:
                continue

            with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
                wtheta = np.load(f)

            # Use gammas from Zarrouk+21
            wtheta_rescaled_zarrouk = wtheta * np.power(bins, -1*(1 - GAMMAS[region][int(rmag_low)]))

            ax.plot(
                bins,
                wtheta_rescaled_zarrouk,
                linewidth=LINEWIDTH,
                linestyle=LINESTYLE[region],
                marker=MARKERS[region],
                color=COLORS[rmag_low],
                #label=f"{rmag_low} < r < {rmag_high}, {region}"
            )

            """
            ax.plot(
                ZARROUK_DATA[int(rmag_low)]["x"],
                ZARROUK_DATA[int(rmag_low)]["y"],
                c="black",
                marker="o",
                linewidth=LINEWIDTH,
                zorder=100, # Large enough to always end up on top
            )
            """

    ax.set_xlim([1e-2, 20])
    ax.set_ylim([1e-4, 10])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    #ax.set_title(f"BGS {brightness} Targets\nComparison with Zarrouk+21")

    ax.set_xscale("log")
    ax.set_yscale("log")

    #ax.legend(loc="upper right")

    ax.grid()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{brightness}_zarrouk_comparison.pdf")

    if SHOW:
        plt.show()
    else:
        plt.clf()


# Plot ratio of each region to the fullsky, for each r-band range
fig, axs = plt.subplots(nrows=len(RMAG_BINS), ncols=1, sharex=True, figsize=(6, 15))
fig.subplots_adjust(hspace=0.21)

fullsky_files = filedicts[ALL_REGIONS]

plotted_lines = []

for ax, (rmag_low, rmag_high) in zip(axs, RMAG_BINS):
    tmp_plotted_lines = []

    with open(f"{PATH_TO_2PCF_FILES}/{ALL_REGIONS}_2PCF_bins_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
        fullsky_bins = np.load(f)
    with open(f"{PATH_TO_2PCF_FILES}/{ALL_REGIONS}_2PCF_wtheta_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
        fullsky_wtheta = np.load(f)

    for region in (BASS_MzLS, DECaLS_NGC, DECaLS_SGC):
        with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_bins_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
            region_bins = np.load(f)
        with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
            region_wtheta = np.load(f)


        line, = ax.plot(
            bins,
            region_wtheta / fullsky_wtheta - 1,
            linewidth=LINEWIDTH,
            linestyle=LINESTYLE[region],
            marker=MARKERS[region],
            color=COLORS[rmag_low],
            #label=f"{rmag_low} < r < {rmag_high}, {region}"
        )
        tmp_plotted_lines.append(line)

        ax.plot(
            [1e-5, 100],
            [0, 0],
            c="black",
        )

    plotted_lines.append(tmp_plotted_lines)

    ax.set_xlim([1e-2, 20])
    ax.set_ylim([-0.2, 0.2])

    ax.set_xscale("log")

    ax.grid()


regions_legend = axs[0].legend(
    [i[0] for i in plotted_lines],
    [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in RMAG_BINS],
    loc="upper left",
    bbox_to_anchor=(1.05, 1.0),
)

_ = axs[0].legend(
    plotted_lines[0],
    [BASS_MzLS, DECaLS_NGC, DECaLS_SGC],
    loc="upper left",
    bbox_to_anchor=(1.05, 0),
)

axs[0].add_artist(regions_legend)

axs[-1].set_xlabel(r"$\theta$ [deg]")
axs[2].set_ylabel(r"$w(\theta)_{region}$ / $w(\theta)_{fullsky} - 1$")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{BRIGHT}_ratios_to_fullsky.pdf")

if SHOW:
    plt.show()
else:
    plt.clf()





# Plot only one rmag bin across 3 regions, ratio to fullsky
fig, ax = plt.subplots(1, 1)
rmag_low, rmag_high = 17, 18
with open(f"{PATH_TO_2PCF_FILES}/{ALL_REGIONS}_2PCF_bins_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
    fullsky_bins = np.load(f)
with open(f"{PATH_TO_2PCF_FILES}/{ALL_REGIONS}_2PCF_wtheta_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
    fullsky_wtheta = np.load(f)

for region in (BASS_MzLS, DECaLS_NGC, DECaLS_SGC):
    with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_bins_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
        region_bins = np.load(f)
    with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{BRIGHT}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
        region_wtheta = np.load(f)

    ax.plot(
        bins,
        region_wtheta / fullsky_wtheta - 1,
        linewidth=LINEWIDTH,
        label=f"{region}"
    )

    ax.plot(
        [1e-5, 100],
        [0, 0],
        c="black",
    )

ax.set_xlim([1e-2, 20])
ax.set_ylim([-0.2, 0.2])

ax.set_xscale("log")

ax.grid()

ax.legend()

ax.set_xlabel(r"$\theta$ [deg]")
ax.set_ylabel(r"$w(\theta)_{region}$ / $w(\theta)_{fullsky} - 1$")
ax.set_title(f"{rmag_low:.1f} < r < {rmag_high:.1f}")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{BRIGHT}_ratios_to_fullsky_pascale.pdf")

if SHOW:
    plt.show()
else:
    plt.clf()


exit() # TODO





# overplot Zarrouk and BASS + DECaLS-NGC + DECaLS-SGC, using gammas
print("Compare all-sky with Zarrouk+21")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    fig, ax = plt.subplots(1, 1)

    region = (BASS_MzLS, DECaLS_NGC, DECaLS_SGC)
    filelist = filedicts[region]

    # Reset colors, to use same colors for any given rmag range across all regions
    ax.set_prop_cycle(None)

    for filename in filelist:
        with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
            bins = np.load(f)

        rmag_low, rmag_high, is_primed = get_rmags_isprimed_from_filename(filename)
        if is_primed:
            continue

        with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_unprimed.npy", "rb") as f:
            wtheta = np.load(f)

        # Use best-fit gammas
        popt, pcov = fit_power_law(bins, wtheta)
        exponent = popt[0]
        gamma = 1 - exponent

        wtheta_rescaled_bestfit = wtheta * np.power(bins, -1*exponent)

        ax.plot(
            bins,
            wtheta_rescaled_bestfit,
            linewidth=LINEWIDTH,
            linestyle=LINESTYLE[region],
            marker=MARKERS[region],
            label=rf"{rmag_low} < r < {rmag_high}, $(\gamma = {gamma:.2f}$ | Z21: $\gamma = {GAMMAS[region][int(rmag_low)]:.2f})$",
        )

        ax.plot(
            ZARROUK_DATA[int(rmag_low)]["x"],
            ZARROUK_DATA[int(rmag_low)]["y"],
            c="black",
            marker="o",
            linewidth=LINEWIDTH,
            zorder=100, # Large enough to always end up on top
        )

    ax.set_xlim([1e-2, 20])
    ax.set_ylim([1e-4, 10])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    ax.set_title(f"BGS {brightness} Targets, Best-Fit Power Law\nComparison with Zarrouk+21")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(loc="upper right")

    ax.grid()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{brightness}_zarrouk_bestfit_comparison.pdf")

    if SHOW:
        plt.show()
    else:
        plt.clf()
