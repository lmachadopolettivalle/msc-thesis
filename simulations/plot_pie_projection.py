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

from sham_model_constants import BLUE, RED

RNG = default_rng(seed=42)

plt.rcParams["font.size"] = "12"
plt.rcParams["figure.figsize"] = (9, 6)

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord
BANDS = ["mag_g", "mag_r", "mag_z"]

DESI_region = directories.DECaLS_NGC

# Range of r apparent magnitude used to select objects
RMAG_MIN, RMAG_MAX = 19, 19.5

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

# TODO determine run_id via some better way, or loop through all existing run_id values
run_id = 146

# Range of theta to be displayed
BEGIN, END = 105, 275

# Load x, y, z positions
# Path to output data from SHAM
SHAM_OUTPUT_PATH = directories.path_interpolation(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_region,
    run_id=run_id,
)

galaxies = {}
for coord in ("x_coord", "y_coord", "z_coord"):
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_{coord}.npy"

    galaxies[coord] = np.load(filename)

# Load magnitudes
for band in BANDS:
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_app_{band}.npy"

    galaxies[band] = np.load(filename)

filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_z.npy"
galaxies["redshift"] = np.load(filename)

# Load whether galaxies are blue or red
filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_blue_red.npy"

galaxies["blue_red"] = np.load(filename)

# Load absolute magnitudes
filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_abs_mag.npy"

galaxies["abs_mag"] = np.load(filename)

# Convert 3D positions into RA, Dec
radii = np.sqrt(galaxies["x_coord"]**2 + galaxies["y_coord"]**2 + galaxies["z_coord"]**2)

theta = np.arccos(galaxies["z_coord"] / radii)
phi = np.arctan2(galaxies["y_coord"], galaxies["x_coord"])

# Note that phi is in range [-pi, pi], but for healpy, must be in range [0, 360 degrees]
phi[phi < 0] += 2 * np.pi

galaxies["RA"] = np.degrees(phi)
galaxies["DEC"] = np.degrees(np.pi/2 - theta)

print(min(galaxies["RA"]), max(galaxies["RA"]))
print(min(galaxies["DEC"]), max(galaxies["DEC"]))

# Filter objects based on r-magnitude
rmag_mask = (galaxies["mag_r"] >= RMAG_MIN) & (galaxies["mag_r"] < RMAG_MAX)
for k, v in galaxies.items():
    galaxies[k] = v[rmag_mask]

# Filter to only galaxies close to declination = 0
mask = (np.abs(galaxies["DEC"]) < 5)
for k, v in galaxies.items():
    galaxies[k] = v[mask]

fraction = 1

ids = RNG.choice(
    len(galaxies["redshift"]),
    size=int(fraction * len(galaxies["redshift"])),
    replace=False,
)

# Create polar plot with galaxies
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

ax.scatter(
    np.radians(galaxies["RA"][ids]),
    galaxies["redshift"][ids],
    s=0.5,
    alpha=0.05,
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
    np.radians(galaxies["RA"][ids]),
    galaxies["redshift"][ids],
    s=0.5,
    alpha=0.1,
    c=colors[ids],
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
# Color galaxies by absolute magnitude
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})

print(min(galaxies["abs_mag"]), max(galaxies["abs_mag"]))

cax = ax.scatter(
    np.radians(galaxies["RA"][ids]),
    galaxies["redshift"][ids],
    s=0.1,
    alpha=1,
    c=galaxies["abs_mag"][ids],
    vmin=-22,
    vmax=-18,
    cmap="gist_rainbow",
)

plt.colorbar(cax, label="Abs. Magnitude")

ax.set_theta_offset(
    np.radians(
        -BEGIN + 90 - (END - BEGIN) / 2
    )
)
ax.set_thetamin(BEGIN)
ax.set_thetamax(END)
ax.set_rmax(0.5)

ax.set_title(f"Galaxy Absolute Magnitudes, {RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}, {DESI_region}")
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

plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/SHAM_polarplot_galaxyabsmag_{RMAG_MIN:.1f}_{RMAG_MAX:.1f}_{DESI_region}.png")

plt.show()
