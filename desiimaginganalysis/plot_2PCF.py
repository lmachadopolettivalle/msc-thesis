from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

from constants import BASS_MzLS, DECaLS_NGC, DECaLS_SGC

# Zarrouk+21 Figure 11, for comparison
SHOW_ZARROUK = False

#REGION = BASS_MzLS
#REGION = DECaLS_NGC
#REGION = DECaLS_SGC
REGION = (BASS_MzLS, DECaLS_NGC, DECaLS_SGC, )

# Use primed or unprimed r-mag data
MAG_R_PRIMED = True
if MAG_R_PRIMED:
    primed_filename_extension = "_primed"
else:
    primed_filename_extension = "_unprimed"

# Data obtained from Figure 11 using https://www.graphreader.com/
ZARROUK_DATA = {
    15: {"x":[0.002,0.005,0.008,0.012,0.018,0.027,0.04,0.057,0.083,0.123,0.177,0.264,0.387,0.563,0.825,1.209,1.785,2.615,3.831,5.613,8.104,12.049,17.523],"y":[0.072,0.09,0.088,0.104,0.106,0.1,0.117,0.114,0.117,0.117,0.116,0.114,0.115,0.104,0.114,0.104,0.11,0.111,0.096,0.07,0.054,0.045,0.023]},
    16: {"x":[0.002,0.005,0.008,0.012,0.018,0.027,0.039,0.058,0.084,0.122,0.179,0.268,0.385,0.568,0.838,1.209,1.785,2.577,3.831,5.451,7.869,11.96,17.268],"y":[0.047,0.051,0.049,0.049,0.048,0.05,0.053,0.053,0.053,0.057,0.058,0.054,0.056,0.056,0.059,0.06,0.065,0.064,0.053,0.042,0.046,0.039,0.027]},
    17: {"x":[0.002,0.005,0.008,0.012,0.018,0.027,0.039,0.058,0.083,0.122,0.182,0.262,0.396,0.576,0.838,1.209,1.798,2.634,3.748,5.491,8.044],"y":[0.028,0.025,0.025,0.026,0.026,0.026,0.026,0.027,0.027,0.027,0.027,0.028,0.027,0.028,0.028,0.027,0.026,0.024,0.018,0.011,0.004]},
    18: {"x":[0.002,0.005,0.008,0.012,0.018,0.026,0.039,0.057,0.084,0.123,0.184,0.266,0.39,0.576,0.831,1.218,1.785,2.577,3.775,5.572,7.585],"y":[0.013,0.013,0.013,0.013,0.014,0.014,0.014,0.015,0.014,0.015,0.014,0.015,0.015,0.015,0.015,0.014,0.013,0.01,0.006,0.004,0.002]},
    19: {"x":[0.002,0.005,0.008,0.012,0.018,0.027,0.039,0.058,0.084,0.123,0.179,0.266,0.387,0.563,0.825,1.218,1.772,2.615,3.803,5.531,8.163],"y":[0.006,0.007,0.007,0.007,0.007,0.008,0.008,0.008,0.008,0.008,0.008,0.009,0.009,0.009,0.009,0.008,0.006,0.005,0.003,0.001,0.001]},
}

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
    (BASS_MzLS, DECaLS_NGC, DECaLS_SGC, ): {
        # NOTE: the paper does not have gammas per r-range for the entire target set,
        # so I am using the DECaLS_SGC gammas, since it is the region with the most targets
        # Could also consider using an average of the gammas for each of the 3 regions instead.
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
FILENAMES = [f for f in FILENAMES if f.startswith(str(REGION)) and ("_bins_" in f) and (f.endswith(f"{primed_filename_extension}.npy"))]

TYPE_FILENAMES = {
    BRIGHT: [f for f in FILENAMES if BRIGHT in f],
    FAINT: [f for f in FILENAMES if FAINT in f],
}

def get_rmag_range_from_filename(filename):
    rmag_range = filename[
        filename.find("range")+len("range"):
        filename.find(f"{primed_filename_extension}.npy")
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

        with open(f"{PATH_TO_2PCF_FILES}/{REGION}_2PCF_wtheta_{type_targets}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}{primed_filename_extension}.npy", "rb") as f:
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

        # Compare against Zarrouk+21 Figure 11
        if SHOW_ZARROUK:
            ax3.scatter(ZARROUK_DATA[int(rmag_low)]["x"], ZARROUK_DATA[int(rmag_low)]["y"])

    for ax in {ax1, ax2, ax3}:
        ax.set_xlim([2e-3, 20])
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(loc="upper right", ncol=2)
        ax.grid()

    ax1.set_title(f"BGS {type_targets} Targets, {REGION}, {'primed' if MAG_R_PRIMED else 'unprimed'}")
    ax2.set_title(f"BGS {type_targets} Targets, {REGION}, {'primed' if MAG_R_PRIMED else 'unprimed'} - Best Fit")
    ax3.set_title(f"BGS {type_targets} Targets, {REGION}, {'primed' if MAG_R_PRIMED else 'unprimed'}\nComparison against Zarrouk+21")

    ax1.set_ylabel(r"$w(\theta)$")
    ax2.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    ax3.set_ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")

    ax1.set_ylim([1e-4, 20])
    ax2.set_ylim([1e-4, 2])
    ax3.set_ylim([1e-4, 2])

    fig1.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}{primed_filename_extension}.pdf")
    #fig2.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}_bestfit{primed_filename_extension}.pdf")
    fig3.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}_zarrouk21{primed_filename_extension}.pdf")

    plt.show()
