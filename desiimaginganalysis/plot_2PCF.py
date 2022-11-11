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
SHOW = False

# Parameters for plotting
plt.rcParams["font.size"] = "12"
LINEWIDTH = 2

BRIGHT = "Bright"
FAINT = "Faint"

LINESTYLE = {
    BASS_MzLS: "solid",
    DECaLS_NGC: "dashed",
    DECaLS_SGC: "dashdot",
    (BASS_MzLS, DECaLS_NGC, DECaLS_SGC): "solid",
}

MARKERS = {
    BASS_MzLS: None,
    DECaLS_NGC: "o",
    DECaLS_SGC: "^",
    (BASS_MzLS, DECaLS_NGC, DECaLS_SGC): "*",
}

COLORS = {
    15: "C0",
    16: "C1",
    17: "C2",
    18: "C3",
    19: "C4",
}

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
        (BASS_MzLS, DECaLS_NGC, DECaLS_SGC, )
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
    (BASS_MzLS, DECaLS_NGC, DECaLS_SGC, ): {
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



# r-primed vs unprimed for BASS and for all-sky
print("r-primed vs unprimed")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    for region in [
        BASS_MzLS,
        (BASS_MzLS, DECaLS_NGC, DECaLS_SGC)
    ]:
        filelist = filedicts[region]

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

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

        ax.set_xlim([2e-3, 20])
        ax.set_ylim([1e-4, 100])
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


# Overplotted curves for the 3 regions, without gamma factor
print("Compare 3 regions without gamma")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    for region in [
        BASS_MzLS,
        DECaLS_NGC,
        DECaLS_SGC,
    ]:
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

            ax.plot(
                bins,
                wtheta,
                linewidth=LINEWIDTH,
                linestyle=LINESTYLE[region],
                marker=MARKERS[region],
                label=f"{rmag_low} < r < {rmag_high}, {region}"
            )

    ax.set_xlim([2e-3, 20])
    ax.set_ylim([1e-4, 100])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$")
    ax.set_title(f"BGS {brightness} Targets\nComparison between Regions (unprimed)")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(loc="upper right", ncol=3)

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

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    for region in [
        BASS_MzLS,
        DECaLS_NGC,
        DECaLS_SGC,
    ]:
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

            # Use gammas from Zarrouk+21
            wtheta_rescaled_zarrouk = wtheta * np.power(bins, -1*(1 - GAMMAS[region][int(rmag_low)]))

            ax.plot(
                bins,
                wtheta_rescaled_zarrouk,
                linewidth=LINEWIDTH,
                linestyle=LINESTYLE[region],
                marker=MARKERS[region],
                label=f"{rmag_low} < r < {rmag_high}, {region}"
            )

            ax.plot(
                ZARROUK_DATA[int(rmag_low)]["x"],
                ZARROUK_DATA[int(rmag_low)]["y"],
                c="black",
                marker="o",
                linewidth=LINEWIDTH,
                zorder=100, # Large enough to always end up on top
            )

    ax.set_xlim([2e-3, 20])
    ax.set_ylim([1e-4, 5])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    ax.set_title(f"BGS {brightness} Targets\nComparison with Zarrouk+21")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(loc="upper right", ncol=3)

    ax.grid()

    #fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{brightness}_zarrouk_comparison.pdf")

    if SHOW:
        plt.show()
    else:
        plt.clf()


# overplot Zarrouk and BASS + DECaLS-NGC + DECaLS-SGC, using gammas
print("Compare all-sky with Zarrouk+21")
for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

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

    ax.set_xlim([2e-3, 20])
    ax.set_ylim([1e-4, 5])
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
