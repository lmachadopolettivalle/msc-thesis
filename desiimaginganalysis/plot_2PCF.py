from matplotlib import pyplot as plt
import numpy as np
import os

REGION = "north"

# Parameters for plotting
plt.rcParams['font.size'] = '12'
LINEWIDTH = 2

# Obtain list of 2PCF files available
PATH_TO_2PCF_FILES = "/cluster/scratch/lmachado/DataProducts/2PCF/"

FILENAMES = sorted(os.listdir(PATH_TO_2PCF_FILES))

# Focus on desired region, on _bins_ files, and on bright vs. faint
FILENAMES = [f for f in FILENAMES if f.startswith(REGION) and "_bins_" in f]

TYPE_FILENAMES = {
    "Bright": [f for f in FILENAMES if "Bright" in f],
    "Faint": [f for f in FILENAMES if "Faint" in f],
}

def get_rmag_range_from_filename(filename):
    rmag_range = filename[
        filename.find("range")+len("range"):
        filename.find(".npy")
    ]
    rmag_low, rmag_high = rmag_range.split("-")
    return float(rmag_low), float(rmag_high)

for type_targets, filelist in TYPE_FILENAMES.items():
    for filename in filelist:
        with open(f"{PATH_TO_2PCF_FILES}/{filename}", "rb") as f:
            bins = np.load(f)

        rmag_low, rmag_high = get_rmag_range_from_filename(filename)

        with open(f"{PATH_TO_2PCF_FILES}/north_2PCF_wtheta_{type_targets}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}.npy", "rb") as f:
            wtheta = np.load(f)

        plt.plot(bins, wtheta, linewidth=LINEWIDTH, label=f"{rmag_low} < r < {rmag_high}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\theta$ [deg]")
    plt.ylabel(r"$w(\theta)$")
    plt.title(f"BGS {type_targets} Targets - {REGION}")
    plt.legend(loc="upper right")

    plt.grid()

    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_2PCF_{type_targets}.pdf")

    plt.show()
