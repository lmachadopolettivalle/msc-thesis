# Using x, y, z coordinates from galaxies
# as outputted by SHAM,
# determine their RA and Dec,
# and display their clustering using Gaussian KDE.

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

import directories

from load_sham_galaxies import load_sham_galaxies

from sham_model_constants import BLUE, RED

# Colors
blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord.
BANDS = ["mag_g", "mag_r", "mag_z"]

# Parameters for plotting
plt.rcParams["figure.figsize"] = (12, 9)
plt.rcParams["font.size"] = "14"
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"

# DESI Region to be loaded
DESI_region = directories.DECaLS_NGC

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"
run_id = 146

# Load galaxy data from SHAM
galaxies = load_sham_galaxies(
    particle_count=PARTICLE_COUNT_PINOCCHIO,
    z_depth=Z_DEPTH,
    pinocchio_region=PINOCCHIO_REGION,
    DESI_region=DESI_region,
    run_id=run_id,
)

####################
# Boundaries of RA and DEC, used to zoom in for visualization
RA_MIN, RA_MAX = 150, 155
DEC_MIN, DEC_MAX = -2.5, 2.5

# Shared Filtering by RA and Dec
ra_dec_filtering = (galaxies["RA"] >= RA_MIN) & (galaxies["RA"] < RA_MAX) & (galaxies["DEC"] >= DEC_MIN) & (galaxies["DEC"] < DEC_MAX)

####################
# Filtering 1 ("Redshift Filtering"):
# - Narrow range in redshift
# - No cut in apparent magnitude
# - Color and size based on absolute magnitude
# - Split into reds vs. blues
REDSHIFT_MIN, REDSHIFT_MAX = 0.3, 0.31

red_mask = (galaxies["blue_red"] == RED)
blue_mask = (galaxies["blue_red"] == BLUE)

redshift_filtering = ra_dec_filtering & (galaxies["redshift"] >= REDSHIFT_MIN) & (galaxies["redshift"] < REDSHIFT_MAX)

redshift_filtering_size = (-240 * galaxies["abs_mag"] - 5020)

blue_redshift_filtering_config = {
    "mask": blue_mask & redshift_filtering,
    "colorbar_label": "Absolute Magnitude",
    "colorbar_field": galaxies["abs_mag"],
    "colormap": "cividis_r",
    "colorbar_vmin": -23,
    "colorbar_vmax": -19,
    "size_field": redshift_filtering_size,
    "title": f"Blue galaxies, {REDSHIFT_MIN:.2f} < z < {REDSHIFT_MAX:.2f}, {DESI_region}",
    "galaxy_color": "blue",
}
red_redshift_filtering_config = {
    "mask": red_mask & redshift_filtering,
    "colorbar_label": "Absolute Magnitude",
    "colorbar_field": galaxies["abs_mag"],
    "colormap": "cividis_r",
    "colorbar_vmin": -23,
    "colorbar_vmax": -19,
    "size_field": redshift_filtering_size,
    "title": f"Red galaxies, {REDSHIFT_MIN:.2f} < z < {REDSHIFT_MAX:.2f}, {DESI_region}",
    "galaxy_color": "red",
}

"""
####################
# Filtering 2 ("App Mag Filtering"):
# - No cut in redshift
# - Cut in apparent magnitude (focus on bright objects)
# - Color and size based on either g-r or r
# - No use of absolute magnitude
# - No split between reds vs. blues
RMAG_MAX = 17

app_mag_filtering = ra_dec_filtering & (galaxies["mag_r"] < RMAG_MAX)

app_mag_filtering_config = {
    "mask": app_mag_filtering,
    "colorbar_label": "g - r",
    "colorbar_field": galaxies["mag_g"] - galaxies["mag_r"],
    "colormap": "RdBu_r",
    "colorbar_vmin": 0.4,
    "colorbar_vmax": 1.0,
    "size_field": np.full(len(galaxies["mag_r"]), 60),
    "title": f"All galaxies, r < {RMAG_MAX:.1f}, {DESI_region} region",
    "galaxy_color": "all",
}
"""

####################
# For each filtering, create following plots:
# - Scatter Plot
# - Gaussian KDE
fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
fig.subplots_adjust(wspace=0.15)
for filtering_config, ax in zip(
        [
            blue_redshift_filtering_config,
            red_redshift_filtering_config,
        ],
        axs
    ):
    # Filter galaxies
    filtered_galaxies = {}
    for k, v in galaxies.items():
        filtered_galaxies[k] = v[filtering_config["mask"]]

    print("Number of objects:", len(filtered_galaxies["mag_r"]))

    scatter = ax.scatter(
        filtered_galaxies["RA"],
        filtered_galaxies["DEC"],
        c=filtering_config["colorbar_field"][filtering_config["mask"]],
        s=filtering_config["size_field"][filtering_config["mask"]],
        vmin=filtering_config["colorbar_vmin"],
        vmax=filtering_config["colorbar_vmax"],
        cmap=filtering_config["colormap"],
    )

    ax.text(
        151.1,
        2.1,
        f"{filtering_config['galaxy_color'].capitalize()} galaxies",
        color=(blue if filtering_config["galaxy_color"] == "blue" else red),
        bbox=dict(facecolor="none", edgecolor="black", pad=6),
        ha="center",
        va="center",
    )

    ax.set_xlabel("RA (degrees)")
    #ax.set_title(filtering_config["title"])

    ax.set_xlim([RA_MIN, RA_MAX])
    ax.set_ylim([DEC_MIN, DEC_MAX])

    ax.set_aspect("equal")

cbar = fig.colorbar(
    scatter,
    label=filtering_config["colorbar_label"],
    ax=axs,
    shrink=0.7,
    fraction=0.046,
    pad=0.04,
)
cbar.ax.invert_yaxis()

axs[0].set_ylabel("Dec (degrees)")

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/Scatter_{DESI_region}.pdf")

plt.show()


####################
# Plot Gaussian KDE
fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
fig.subplots_adjust(wspace=0.15)
for filtering_config, ax in zip(
        [
            blue_redshift_filtering_config,
            red_redshift_filtering_config,
        ],
        axs
    ):

    # Filter galaxies
    filtered_galaxies = {}
    for k, v in galaxies.items():
        filtered_galaxies[k] = v[filtering_config["mask"]]

    print("Number of objects:", len(filtered_galaxies["mag_r"]))

    X, Y = np.mgrid[RA_MIN:RA_MAX:128j, DEC_MIN:DEC_MAX:128j]

    positions = np.vstack([X.ravel(), Y.ravel()])

    values = np.vstack([
        filtered_galaxies["RA"],
        filtered_galaxies["DEC"],
    ])

    kernel = stats.gaussian_kde(
        values,
        bw_method=0.1,
    )

    Z = np.reshape(kernel(positions).T, X.shape)

    cax = ax.imshow(
        np.rot90(Z),
        cmap="jet",
        extent=[RA_MIN, RA_MAX, DEC_MIN, DEC_MAX],
        vmin=0,
        vmax=0.6,
    )

    ax.set_xlim([RA_MIN, RA_MAX])
    ax.set_ylim([DEC_MIN, DEC_MAX])

    ax.set_xlabel("RA (degrees)")
    #ax.set_title(filtering_config["title"])

    ax.set_aspect("equal")

axs[0].set_ylabel("Dec (degrees)")
fig.colorbar(
    cax,
    label=r"PDF ($deg^{-2}$)",
    ax=axs,
    shrink=0.7,
    fraction=0.046,
    pad=0.04,
)

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/GaussianKDE_{DESI_region}.pdf")

plt.show()
