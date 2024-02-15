# Compare 2PCF with extinction ON and OFF

from matplotlib import pyplot as plt
import numpy as np
import os

from constants import *

# Parameters for plotting
plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"
LINEWIDTH = 2

BRIGHT = "Bright"
FAINT = "Faint"

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
        DECaLS_SGC,
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
    for brightness in [BRIGHT]
}

region = DECaLS_SGC

for brightness, filedicts in TYPE_FILENAMES.items():
    # Focus on BRIGHT
    if brightness != BRIGHT:
        continue

    plotted_lines = []

    fig, ax = plt.subplots(1, 1)

    filelist = filedicts[region]

    # Reset colors, to use same colors for any given rmag range across all regions
    ax.set_prop_cycle(None)

    for rmag_low, rmag_high in RMAG_BINS:
        tmp_plotted_lines = []
        for is_extinction in [True, False]:
            # Only computed noextinction for primed, and extinction for unprimed
            # But for DECaLS primed and unprimed are the same
            is_primed = not is_extinction
            with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_bins_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if is_primed else 'unprimed'}{'_noextinction' if not is_extinction else ''}.npy", "rb") as f:
                bins = np.load(f)
            with open(f"{PATH_TO_2PCF_FILES}/{region}_2PCF_wtheta_{brightness}_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_{'primed' if is_primed else 'unprimed'}{'_noextinction' if not is_extinction else ''}.npy", "rb") as f:
                wtheta = np.load(f)

            line, = ax.plot(
                bins,
                wtheta,
                linewidth=LINEWIDTH,
                color=COLORS[rmag_low],
                linestyle=("solid" if is_extinction else "dashed"),
                #label=f"{rmag_low} < r < {rmag_high}",
            )
            tmp_plotted_lines.append(line)
        plotted_lines.append(tmp_plotted_lines)

    ax.set_xlim([1e-1, 1])
    ax.set_ylim([1e-3, 10])
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$w(\theta)$")
    #ax.set_title(f"BGS {brightness} Targets\nComparison between Regions (unprimed)")

    ax.set_xscale("log")
    ax.set_yscale("log")

    rmags_legend = ax.legend(
        [i[0] for i in plotted_lines],
        [f"{rmag_low:.1f} < r < {rmag_high:.1f}" for (rmag_low, rmag_high) in RMAG_BINS],
        loc="upper right",
    )

    _ = ax.legend(
        plotted_lines[0],
        ["Extinction", "No Extinction"],
        loc="lower left",
    )

    ax.add_artist(rmags_legend)

    ax.grid()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/2PCF_{region}_noextinction_comparison.pdf")

    plt.show()
