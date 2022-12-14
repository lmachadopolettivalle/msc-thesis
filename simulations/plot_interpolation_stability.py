# For a given parameter space run,
# plot absolute magnitude vs. halo/subhalo mass
# for several different redshifts.
# Make two plots, one for blue and one for red galaxies.

from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from tqdm import tqdm

from manage_parameter_space import get_details_of_run

plt.rcParams["font.size"] = "12"

# TODO find better way to obtain the run_id,
# or loop through all available run IDs
run_id = 100
pinocchio_particle_count = 2048

infile_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{pinocchio_particle_count}cubed/{run_id}/interpolation_outputs/"

infile_hist_blue = "blue_lim_interp.npz"
infile_hist_red = "red_lim_interp.npz"

run_details = get_details_of_run(run_id)
# SHAM parameters
M_limit = run_details["mass_cut"] # mass limit for assigning blue or red galaxies to halos, [Msun/h]

# ----------------------------------------------------–
# LOAD HALO-SUBHALO HISTOGRAMS
# -----------------------------------------------------
print("Loading blue galaxies")
with np.load(infile_dir + infile_hist_blue) as data:
    lim_abs_mag_blue = data["lim_abs_mag_blue"]
    blue_bin_edges_z = data["bin_edges_z"]
    blue_bin_edges_mass = data["bin_edges_mass"]

print("Loading red galaxies")
with np.load(infile_dir + infile_hist_red) as data:
    lim_abs_mag_red = data["lim_abs_mag_red"]
    red_bin_edges_z = data["bin_edges_z"]
    red_bin_edges_mass = data["bin_edges_mass"]

# Create arrays for interpolation
blue_mass_edges_stacked = np.reshape(
    np.tile(blue_bin_edges_mass, len(blue_bin_edges_z)),
    (len(blue_bin_edges_mass) * len(blue_bin_edges_z), ),
)
red_mass_edges_stacked = np.reshape(
    np.tile(red_bin_edges_mass, len(red_bin_edges_z)),
    (len(red_bin_edges_mass) * len(red_bin_edges_z), ),
)

blue_z_edges_stacked = np.reshape(
    np.ndarray.flatten(
        np.transpose(
            np.tile(
                blue_bin_edges_z,
                (len(blue_bin_edges_mass), 1),
            )
        )
    ),
    (len(blue_bin_edges_mass)*len(blue_bin_edges_z), ),
)
red_z_edges_stacked = np.reshape(
    np.ndarray.flatten(
        np.transpose(
            np.tile(
                red_bin_edges_z,
                (len(red_bin_edges_mass), 1),
            )
        )
    ),
    (len(red_bin_edges_mass)*len(red_bin_edges_z), ),
)

# Choose input values for mass and redshift
input_masses = np.linspace(11, 15, num=100)
input_redshifts = np.linspace(0, 0.5, 10)

fig_blue, ax_blue = plt.subplots(1, 1, figsize=(8, 6))
fig_red, ax_red = plt.subplots(1, 1, figsize=(8, 6))

# Interpolate for each redshift
for z in tqdm(input_redshifts):
    blue_absmags = griddata(
        (blue_z_edges_stacked, np.log10(blue_mass_edges_stacked)),
        np.ndarray.flatten(lim_abs_mag_blue),
        (z, input_masses),
    )
    red_absmags = griddata(
        (red_z_edges_stacked, np.log10(red_mass_edges_stacked)),
        np.ndarray.flatten(lim_abs_mag_red),
        (z, input_masses),
    )

    ax_blue.plot(
        10**input_masses,
        blue_absmags,
        label=f"z = {z:.2f}",
        linewidth=2,
    )
    ax_red.plot(
        10**input_masses,
        red_absmags,
        label=f"z = {z:.2f}",
        linewidth=2,
    )

ax_blue.set_title(f"Blue, mass cut = {M_limit:.1e}")
ax_red.set_title(f"Red, mass cut = {M_limit:.1e}")

for ax in {ax_blue, ax_red}:
    ax.set_xlabel(r"Halo/Subhalo Mass ($M_{\odot}/h$)")
    ax.set_ylabel("Galaxy's Absolute Magnitude")

    ax.invert_yaxis()

    ax.set_xscale("log")

    ax.set_ylim([-17, -23])

    ax.legend(loc="lower right")

    ax.grid()

fig_blue.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/blue_interpolation.pdf")
fig_red.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/red_interpolation.pdf")

plt.show()
