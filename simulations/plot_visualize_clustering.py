# Using x, y, z coordinates from galaxies
# as outputted by SHAM,
# determine their RA and Dec,
# and display their clustering using Gaussian KDE.

import mpl_scatter_density # Following guide from https://stackoverflow.com/a/64105308
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

import directories

from sham_model_constants import BLUE, RED

# Parameters for plotting
plt.rcParams["font.size"] = "12"
plt.rcParams["figure.figsize"] = (8, 6)

# DESI Region to be loaded
DESI_region = directories.BASS_MzLS

# Boundaries of RA and DEC, used to zoom in for visualization
RA_MIN, RA_MAX = 150, 180
DEC_MIN, DEC_MAX = 40, 60

# Range of r-magnitude used to filter galaxies
RMAG_MIN, RMAG_MAX = 14, 19.5

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord.
BANDS = ["mag_g", "mag_r", "mag_z"]

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

# TODO determine run_id via some better way, or loop through all existing run_id values
run_id = 146

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


# Load whether galaxies are blue or red
filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_blue_red.npy"

galaxies["blue_red"] = np.load(filename)

# Convert 3D positions into RA, Dec
radii = np.sqrt(galaxies["x_coord"]**2 + galaxies["y_coord"]**2 + galaxies["z_coord"]**2)

theta = np.arccos(galaxies["z_coord"] / radii)
phi = np.arctan2(galaxies["y_coord"], galaxies["x_coord"])

# Note that phi is in range [-pi, pi], but for healpy, must be in range [0, 360 degrees]
phi[phi < 0] += 2 * np.pi

galaxies["RA"] = np.degrees(phi)
galaxies["DEC"] = np.degrees(np.pi/2 - theta)

# Filter galaxies based on RA and DEC boundaries
ra_dec_filter = (galaxies["RA"] >= RA_MIN) & (galaxies["RA"] < RA_MAX) & (galaxies["DEC"] >= DEC_MIN) & (galaxies["DEC"] < DEC_MAX)
for k, v in galaxies.items():
    galaxies[k] = v[ra_dec_filter]

# Filter galaxies based on r-magnitude
rmag_filter = (galaxies["mag_r"] >= RMAG_MIN) & (galaxies["mag_r"] < RMAG_MAX)
for k, v in galaxies.items():
    galaxies[k] = v[rmag_filter]

# Separate red and blue galaxies
red_mask = galaxies["blue_red"] == RED
blue_mask = galaxies["blue_red"] == BLUE
red_galaxies = {k: v[red_mask] for k, v in galaxies.items()}
blue_galaxies = {k: v[blue_mask] for k, v in galaxies.items()}

# Plot 2D histogram using mpl scatter
for g, color in (
        (blue_galaxies, "Blue"),
        (red_galaxies, "Red"),
    ):

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection":"scatter_density"})

    density = ax.scatter_density(
        g["RA"],
        g["DEC"],
        vmin=0,
        vmax=10,
        #vmax=5,
        cmap=plt.cm.Greys,
    )
    fig.colorbar(
        density,
        label="Number of points per pixel",
    )

    ax.set_xlabel("RA (degrees)")
    ax.set_ylabel("Dec (degrees)")
    ax.set_title(f"{color} galaxies, {DESI_region} region\n{RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}")

    ax.set_xlim([RA_MIN, RA_MAX])
    ax.set_ylim([DEC_MIN, DEC_MAX])

    plt.tight_layout()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/ScatterDensity_{DESI_region}_{RMAG_MIN}_{RMAG_MAX}_{color}.pdf")

    plt.show()

# Plot Gaussian KDE
SIZE = 35000
X, Y = np.mgrid[RA_MIN:RA_MAX:100j, DEC_MIN:DEC_MAX:100j]

positions = np.vstack([X.ravel(), Y.ravel()])

for g, color in (
        (blue_galaxies, "Blue"),
        (red_galaxies, "Red"),
    ):
    values = np.vstack([
        g["RA"][:SIZE],
        g["DEC"][:SIZE],
    ])

    kernel = stats.gaussian_kde(values)

    Z = np.reshape(kernel(positions).T, X.shape)

    fig, ax = plt.subplots()

    cax = ax.imshow(
        np.rot90(Z),
        cmap=plt.cm.Greys,
        extent=[RA_MIN, RA_MAX, DEC_MIN, DEC_MAX],
        vmin=0,
        vmax=5e-3,
        #vmax=3e-3,
    )

    fig.colorbar(cax)

    ax.set_xlim([RA_MIN, RA_MAX])
    ax.set_ylim([DEC_MIN, DEC_MAX])

    ax.set_xlabel("RA (degrees)")
    ax.set_ylabel("Dec (degrees)")
    ax.set_title(f"Gaussian KDE\n{color} galaxies, {DESI_region} region\n{RMAG_MIN:.1f} < r < {RMAG_MAX:.1f}")

    plt.tight_layout()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/GaussianKDE_{DESI_region}_{RMAG_MIN}_{RMAG_MAX}_{color}.pdf")

    plt.show()
