# For a given parameter space run,
# plot absolute magnitude vs. halo/subhalo mass
# for several different redshifts.
# Make two plots, one for blue and one for red galaxies.

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import re
from tqdm import tqdm

from manage_parameter_space import get_details_of_run

plt.rcParams["font.size"] = "12"

pinocchio_particle_count = 2048

DESIRED_RUN_IDS = [100] + list(range(106, 139 + 1))

# Store number of galaxies generated in each run
galaxy_counts = {} # Key = (num_z_bins, num_mass_bins), value = {"total": 0, "blue": 0, "red": 0}
z_bins = set()
mass_bins = set()

def extract_galaxy_counts(filename):
    with open(filename, 'r') as f:
        text = f.read()

    total = int(
        re.search("total number of galaxies after magnitude cuts = (.+)\n", text).groups()[0]
    )
    blue = int(float(
        re.search("number of blue galaxies after magnitude cuts = (.+)\n", text).groups()[0]
    ))
    red = int(float(
        re.search("number of red galaxies after magnitude cuts = (.+)\n", text).groups()[0]
    ))

    return (
        total,
        blue,
        red,
    )

# Loop through desired run IDs
for run_id in DESIRED_RUN_IDS:
    # File with output of SHAM run.
    # Used to determine effective mass limit used in SHAM
    sham_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/sham_int_job_{pinocchio_particle_count}_{run_id}_output"

    # Obtain details of this SHAM run,
    # in particular the desired mass limit
    run_details = get_details_of_run(run_id)
    num_z_bins = run_details["num_z_bins"]
    num_mass_bins = run_details["num_mass_bins"]

    z_bins.add(num_z_bins)
    mass_bins.add(num_mass_bins)

    total, blue, red = extract_galaxy_counts(sham_output_filename)

    galaxy_counts[(num_z_bins, num_mass_bins)] = {
        "total": total,
        "blue": blue,
        "red": red,
    }

z_bins = sorted(z_bins)
mass_bins = sorted(mass_bins)

# Row keeps mass bin constant, varying z bins
# Column keeps z bin constant, varying mass bins
for k in ["total", "blue", "red"]:
    grid = [[galaxy_counts[(z, m)][k] for z in sorted(z_bins)] for m in sorted(mass_bins)]

    fig, ax = plt.subplots(1, 1, figsize=(12, 9))

    im = ax.imshow(grid, cmap=mpl.cm.RdBu_r)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(z_bins)), labels=z_bins)
    ax.set_yticks(np.arange(len(mass_bins)), labels=mass_bins)

    ax.set_xlabel("Mass bins")
    ax.set_ylabel("Redshift bins")

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Create colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Galaxy Count", rotation=-90, va="bottom")

    # Loop over data dimensions and create text annotations.
    for i in range(len(mass_bins)):
        for j in range(len(z_bins)):
            ax.text(j, i, grid[i][j], ha="center", va="center", color="black")

    ax.set_title(f"{k.capitalize()} galaxy counts, different parameters")

    fig.tight_layout()

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/galaxy_counts_{k}_from_interpolation_{pinocchio_particle_count}.pdf")

    plt.show()

    # Close figures to reduce memory usage,
    # since pyplot keeps figures in memory until
    # the end of the program by default
    plt.close(fig)
