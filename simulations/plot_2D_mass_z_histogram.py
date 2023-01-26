import matplotlib
from matplotlib import pyplot as plt
import numpy as np

from manage_parameter_space import get_details_of_run

plt.rcParams["font.size"] = "12"
#plt.rcParams["figure.figsize"] = (30, 8) # Default: (6.4, 4.8)

pinocchio_particle_count = 2048

run_id = 143

# Directory containing output data from SHAM
infile_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{pinocchio_particle_count}cubed/{run_id}/interpolation_outputs/"
infile_hist_blue = "blue_lim_interp.npz"
infile_hist_red = "red_lim_interp.npz"

data = {}

with np.load(infile_dir + infile_hist_blue) as d:
    data["blue"] = {
        "lim_abs_mag": d["lim_abs_mag_blue"],
        "bin_edges_z": d["bin_edges_z"],
        "bin_edges_mass": d["bin_edges_mass"],
    }

with np.load(infile_dir + infile_hist_red) as d:
    data["red"] = {
        "lim_abs_mag": d["lim_abs_mag_red"],
        "bin_edges_z": d["bin_edges_z"],
        "bin_edges_mass": d["bin_edges_mass"],
    }

for color in ("blue", "red"):
    lim_abs_mag = data[color]["lim_abs_mag"]
    bin_edges_z = data[color]["bin_edges_z"]
    bin_edges_mass = data[color]["bin_edges_mass"]

    plt.pcolormesh(bin_edges_z, bin_edges_mass, lim_abs_mag.T, vmin=-24, vmax=-18)

    cbar = plt.colorbar()
    cbar.ax.invert_yaxis()

    plt.yscale("log")
    plt.ylim([5e10, 1e15])

    plt.xlabel("z")
    plt.ylabel(r"Halo/Subhalo mass [M$_{\odot}$/h]")

    plt.title(f"Absolute Magnitudes of {color.capitalize()} Galaxies")

    plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/{color}_run{run_id}_2D_histogram.pdf")

    plt.show()
