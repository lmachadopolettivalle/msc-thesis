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
from numpy.random import default_rng
import os
import pandas as pd
from tqdm import tqdm

import directories

from load_sham_galaxies import load_sham_galaxies

from sham_model_constants import BLUE, RED

RNG = default_rng(seed=420)

plt.rcParams["font.size"] = "12"
plt.rcParams["figure.figsize"] = (18, 12)

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
dec_mask = (np.abs(galaxies["DEC"]) < 0.1)
for k, v in galaxies.items():
    galaxies[k] = v[dec_mask]

####################
# Create polar plot with galaxies
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

ax.scatter(
    np.radians(galaxies["RA"]),
    galaxies["redshift"],
    s=0.5,
    alpha=1,
    c="black",
)

ax.set_theta_offset(
    np.radians(
        -BEGIN + 90 - (END - BEGIN) / 2
    )
)
ax.set_thetamin(BEGIN)
ax.set_thetamax(END)
ax.set_rmax(0.5)

ax.set_title(f"Galaxy Positions, {RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}, {DESI_region}")
ax.text(
    np.radians(90),
    0.3,
    "Redshift",
    rotation=(90 - (END - BEGIN)/2),
    ha="center",
    va="center",
)

ax.text(
    np.radians((BEGIN + END) / 2),
    0.6,
    "RA [degrees]",
    ha="center",
    va="center",
)

fig.subplots_adjust(bottom=0, top=1, left=0, right=1)

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/SHAM_polarplot_galaxypositions_{RMAG_MIN:.1f}_{RMAG_MAX:.1f}_{DESI_region}.png")
plt.show()

####################
# Color galaxies by red vs. blue
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

blue_color = "#1f77b4"
red_color = "#d62728"

colors = np.full(len(galaxies["blue_red"]), red_color)

colors[
    galaxies["blue_red"] == BLUE
] = blue_color

ax.scatter(
    np.radians(galaxies["RA"]),
    galaxies["redshift"],
    s=0.5,
    alpha=0.3,
    c=colors,
)

ax.set_theta_offset(
    np.radians(
        -BEGIN + 90 - (END - BEGIN) / 2
    )
)
ax.set_thetamin(BEGIN)
ax.set_thetamax(END)
ax.set_rmax(0.5)

ax.set_title(f"Galaxy Red vs. Blue, {RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}, {DESI_region}")
ax.text(
    np.radians(90),
    0.3,
    "Redshift",
    rotation=(90 - (END - BEGIN)/2),
    ha="center",
    va="center",
)

ax.text(
    np.radians((BEGIN + END) / 2),
    0.6,
    "RA [degrees]",
    ha="center",
    va="center",
)

fig.subplots_adjust(bottom=0, top=1, left=0, right=1)

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/SHAM_polarplot_galaxyredvsblue_{RMAG_MIN:.1f}_{RMAG_MAX:.1f}_{DESI_region}.png")

plt.show()

####################
# Color galaxies by g - r
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

print(min(galaxies["abs_mag"]), max(galaxies["abs_mag"]))

cax = ax.scatter(
    np.radians(galaxies["RA"]),
    galaxies["redshift"],
    s=0.1,
    alpha=1,
    c=(galaxies["mag_g"] - galaxies["mag_r"]),
    vmin=0.4,
    vmax=1,
    cmap="RdBu_r",
)

plt.colorbar(cax, label="g - r", fraction=0.046, pad=0.04, shrink=0.5)

ax.set_theta_offset(
    np.radians(
        -BEGIN + 90 - (END - BEGIN) / 2
    )
)
ax.set_thetamin(BEGIN)
ax.set_thetamax(END)
ax.set_rmax(0.5)

ax.set_title(f"Galaxy g - r Colors, {RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}, {DESI_region}")
ax.text(
    np.radians(90),
    0.3,
    "Redshift",
    rotation=(90 - (END - BEGIN)/2),
    ha="center",
    va="center",
)

ax.text(
    np.radians((BEGIN + END) / 2),
    0.6,
    "RA [degrees]",
    ha="center",
    va="center",
)

fig.subplots_adjust(bottom=0, top=0.8, left=0.1, right=0.7)

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/SHAM_polarplot_gminusr_{RMAG_MIN:.1f}_{RMAG_MAX:.1f}_{DESI_region}.png")

plt.show()
