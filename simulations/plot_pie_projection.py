# GOALS:
# https://www.astroexplorer.org/details/apj391584f1
# https://www.astroexplorer.org/details/apj391584f2
# https://www.astroexplorer.org/details/apj391584f3


# Using x, y, z coordinates from galaxies
# as outputted by SHAM,
# determine their RA and Dec,
# and display them in a polar projection plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import directories

from load_sham_galaxies import load_sham_galaxies

from sham_model_constants import BLUE, RED

plt.rcParams["font.size"] = "14"
plt.rcParams["figure.figsize"] = (18, 14)

DESI_region = directories.DECaLS_NGC

# Range of r apparent magnitude used to select objects
RMAG_MIN, RMAG_MAX = -np.inf, 19.5

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"
run_id = 146

# Range of theta to be displayed
BEGIN, END = 105, 275

# Load galaxy data from SHAM
galaxies = load_sham_galaxies(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_region,
    run_id=run_id,
)

print(min(galaxies["RA"]), max(galaxies["RA"]))
print(min(galaxies["DEC"]), max(galaxies["DEC"]))

####################
# Separate red and blue galaxies
red_mask = galaxies["blue_red"] == RED
blue_mask = galaxies["blue_red"] == BLUE
red_galaxies = {k: v[red_mask] for k, v in galaxies.items()}
blue_galaxies = {k: v[blue_mask] for k, v in galaxies.items()}

####################
# Filter objects based on r-magnitude
rmag_mask = (galaxies["mag_r"] >= RMAG_MIN) & (galaxies["mag_r"] < RMAG_MAX)
for k, v in galaxies.items():
    galaxies[k] = v[rmag_mask]

# Filter to only galaxies close to declination = 0
dec_mask = (np.abs(galaxies["DEC"]) < 0.15)
for k, v in galaxies.items():
    galaxies[k] = v[dec_mask]

####################
# Configurations for the different pie plots
blue_color = "#1f77b4"
red_color = "#d62728"

bluered_colors = np.full(len(galaxies["blue_red"]), red_color)

bluered_colors[
    galaxies["blue_red"] == BLUE
] = blue_color


simple_scatter_config = {
    "colors": "black",
    "size": 0.75,
    "vmin": None,
    "vmax": None,
    "colormap": None,
    "colorbar_label": None,
    "title": f"Galaxy Positions, r < {RMAG_MAX:.1f}, {DESI_region}",
    "filename": f"SHAM_polarplot_galaxypositions_{RMAG_MAX:.1f}_{DESI_region}.png",
}
bluered_config = {
    "colors": bluered_colors,
    "size": 3,
    "vmin": None,
    "vmax": None,
    "colormap": None,
    "colorbar_label": None,
    "title": f"Galaxy Red vs. Blue, r < {RMAG_MAX:.1f}, {DESI_region}",
    "filename": f"SHAM_polarplot_galaxyredvsblue_{RMAG_MAX:.1f}_{DESI_region}.png",
}
gminusr_config = {
    "colors": (galaxies["mag_g"] - galaxies["mag_r"]),
    "size": 1.25,
    "vmin": 0.4,
    "vmax": 1,
    "colormap": "RdBu_r",
    "colorbar_label": "g - r",
    "title": f"Galaxy g - r Colors, r < {RMAG_MAX:.1f}, {DESI_region}",
    "filename": f"SHAM_polarplot_gminusr_{RMAG_MAX:.1f}_{DESI_region}.png",
}

for config in [
        simple_scatter_config,
        bluered_config,
        gminusr_config,
    ]:
    # Create polar plot with galaxies
    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

    cax = ax.scatter(
        np.radians(galaxies["RA"]),
        galaxies["redshift"],
        s=config["size"],
        c=config["colors"],
        vmin=config["vmin"],
        vmax=config["vmax"],
        cmap=config["colormap"],
    )

    if config["colormap"] is not None:
        plt.colorbar(cax, label=config["colorbar_label"], fraction=0.046, pad=0.04, shrink=0.5)

    ax.set_theta_offset(
        np.radians(
            -BEGIN + 90 - (END - BEGIN) / 2
        )
    )
    ax.set_thetamin(BEGIN)
    ax.set_thetamax(END)
    ax.set_rmax(0.3)

    ax.set_title(config["title"])
    ax.text(
        np.radians(90),
        0.15,
        "Redshift",
        ha="center",
        va="center",
    )

    ax.text(
        np.radians((BEGIN + END) / 2),
        0.35,
        "RA [degrees]",
        ha="center",
        va="center",
    )

    if config["colormap"] is None:
        fig.subplots_adjust(bottom=0, top=1, left=0.1, right=1)
    else:
        fig.subplots_adjust(bottom=0, top=0.9, left=0.1, right=0.8)

    plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/{config['filename']}")
    plt.show()
