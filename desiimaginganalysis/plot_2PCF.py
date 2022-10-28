from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

from constants import BASS_MzLS, DECaLS_NGC, DECaLS_SGC


#REGION = BASS_MzLS
REGION = DECaLS_NGC
#REGION = DECaLS_SGC

BRIGHT = "Bright"
FAINT = "Faint"

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
}

# Parameters for plotting
plt.rcParams['font.size'] = '12'
LINEWIDTH = 2

# Obtain list of 2PCF files available
PATH_TO_2PCF_FILES = "/cluster/scratch/lmachado/DataProducts/2PCF/"

FILENAMES = sorted(os.listdir(PATH_TO_2PCF_FILES))

# Focus on desired region, on _bins_ files, and on bright vs. faint
FILENAMES = [f for f in FILENAMES if f.startswith(REGION) and "_bins_" in f]

TYPE_FILENAMES = {
    BRIGHT: [f for f in FILENAMES if BRIGHT in f],
    FAINT: [f for f in FILENAMES if FAINT in f],
}

def get_rmag_range_from_filename(filename):
    rmag_range = filename[
        filename.find("range")+len("range"):
        filename.find(".npy")
    ]
    rmag_low, rmag_high = rmag_range.split("-")
    return float(rmag_low), float(rmag_high)

# Create figures
# 1. without any gamma factor
# 2. with our best-fit gamma factor
# 3. with gamma from Zarrouk+2021
fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
fig3, ax3 = plt.subplots(1, 1, figsize=(8, 6))
# Change colors used in plot
for ax in {ax1, ax2, ax3}:
    ax.set_prop_cycle("color", ["orange", "purple", "red", "green", "blue"])

for type_targets, filelist in TYPE_FILENAMES.items():
    # Focus on BRIGHT targets
    if type_targets != BRIGHT:
        continue

    for filename in filelist:
        with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
            bins = np.load(f)

        rmag_low, rmag_high = get_rmag_range_from_filename(filename)

        with open(f"{PATH_TO_2PCF_FILES}/{REGION}_2PCF_wtheta_{type_targets}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}.npy", "rb") as f:
            wtheta = np.load(f)

        # Fit 2PCF to power law of theta
        popt, pcov = fit_power_law(bins, wtheta)
        exponent = popt[0]
        gamma = 1 - exponent
        print(f"Exponent for r range {rmag_low}-{rmag_high} is {exponent}, which means gamma = {gamma}")
        print(pcov)

        wtheta_rescaled_bestfit = wtheta * np.power(bins, -1*exponent)
        wtheta_rescaled_zarrouk = wtheta * np.power(bins, -1*(1 - GAMMAS[REGION][int(rmag_low)]))

        ax1.plot(bins, wtheta, linewidth=LINEWIDTH, label=f"{rmag_low} < r < {rmag_high}")
        ax2.plot(bins, wtheta_rescaled_bestfit, linewidth=LINEWIDTH, label=f"{rmag_low} < r < {rmag_high}" + " " + rf"$(\gamma = {gamma:.2f})$")
        ax3.plot(bins, wtheta_rescaled_zarrouk, linewidth=LINEWIDTH, label=f"{rmag_low} < r < {rmag_high}" + " " + rf"$(\gamma = {GAMMAS[REGION][int(rmag_low)]:.2f})$")

    for ax in {ax1, ax2, ax3}:
        ax.set_xlim([2e-3, 20])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.legend(loc="upper right", ncol=2)
        ax.grid()

    ax1.set_title(f"BGS {type_targets} Targets, {REGION}")
    ax2.set_title(f"BGS {type_targets} Targets, {REGION} - Best Fit")
    ax3.set_title(f"BGS {type_targets} Targets, {REGION}\nComparison against Zarrouk+21")

    ax1.set_ylabel(r"$w(\theta)$")
    ax2.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    ax3.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")

    ax1.set_ylim([1e-4, 20])
    ax2.set_ylim([1e-4, 2])
    ax3.set_ylim([1e-4, 2])

    fig1.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}.pdf")
    #fig2.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}_bestfit.pdf")
    fig3.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}_zarrouk21.pdf")

    plt.show()
