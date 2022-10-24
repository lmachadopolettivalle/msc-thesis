from matplotlib import pyplot as plt
import numpy as np
import os

REGION = "north"

BRIGHT = "Bright"
FAINT = "Faint"

# Gammas from Table 2 of Zarrouk+2021
GAMMAS = {
    "north": { # "north" refers to the BASS/MzLS rows
        15: 1.642,
        16: 1.744,
        17: 1.776,
        18: 1.750,
        19: 1.725,
    },
    "south": { # "south" refers to the DECaLS-SGC rows. NOTE I know part of our south targets are in the DECaLS-NGC, but it's a bit annoying to split the south for this gamma comparison.
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

for type_targets, filelist in TYPE_FILENAMES.items():
    # Focus on BRIGHT targets
    if type_targets != BRIGHT:
        continue

    for filename in filelist:
        with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
            bins = np.load(f)

        rmag_low, rmag_high = get_rmag_range_from_filename(filename)

        with open(f"{PATH_TO_2PCF_FILES}/north_2PCF_wtheta_{type_targets}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}.npy", "rb") as f:
            wtheta = np.load(f)

        wtheta_rescaled = wtheta * np.power(bins, -(1-GAMMAS[REGION][int(rmag_low)]))

        plt.plot(bins, wtheta_rescaled, linewidth=LINEWIDTH, label=f"{rmag_low} < r < {rmag_high}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\theta$ [deg]")
    plt.ylabel(r"$w(\theta)$ x $\theta^{-(1 - \gamma)}$")
    plt.title(f"BGS {type_targets} Targets - {REGION}")
    plt.legend(loc="upper right")

    plt.grid()

    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}.pdf")

    plt.show()
