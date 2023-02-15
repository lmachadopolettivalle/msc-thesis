# For a given parameter space run,
# plot absolute magnitude vs. halo/subhalo mass
# for several different redshifts.
# Make two plots, one for blue and one for red galaxies.

from matplotlib import pyplot as plt
import numpy as np
import re
from scipy.interpolate import griddata
from tqdm import tqdm

import directories

from manage_parameter_space import get_details_of_run

plt.rcParams["font.size"] = "12"

pinocchio_particle_count = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

DESI_region = directories.BASS_MzLS

# Loop through desired run IDs
DESIRED_RUN_IDS = [146]

for run_id in tqdm(DESIRED_RUN_IDS):
    # File with output of SHAM run.
    # Used to determine effective mass limit used in SHAM
    sham_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/sham_int_job_{pinocchio_particle_count}_{run_id}_output"

    # Directory containing output data from SHAM
    infile_dir = directories.path_interpolation(
        particle_count=pinocchio_particle_count,
        z_depth=Z_DEPTH,
        pinocchio_region=PINOCCHIO_REGION,
        DESI_region=DESI_region,
        run_id=run_id,
    )
    infile_hist_blue = "blue_lim_interp.npz"
    infile_hist_red = "red_lim_interp.npz"

    # Obtain details of this SHAM run,
    # in particular the desired mass limit
    run_details = get_details_of_run(run_id)
    M_limit = run_details["mass_cut"] # mass limit for assigning blue or red galaxies to halos, [Msun/h]
    num_z_bins = run_details["num_z_bins"]
    num_mass_bins = run_details["num_mass_bins"]

    # Obtain effective mass limit used
    with open(sham_output_filename, 'r') as f:
        text = f.read()

    M_limit_effective = float(
        re.search("M_limit_effective\s+=\s+(.+)\n", text).groups()[0]
    )

    # ----------------------------------------------------–
    # LOAD HALO-SUBHALO HISTOGRAMS
    # -----------------------------------------------------
    with np.load(infile_dir + infile_hist_blue) as data:
        lim_abs_mag_blue = data["lim_abs_mag_blue"]
        blue_bin_edges_z = data["bin_edges_z"]
        blue_bin_edges_mass = data["bin_edges_mass"]

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
    input_redshifts = np.linspace(0.05, 0.5, 5)

    fig_blue, ax_blue = plt.subplots(1, 1, figsize=(8, 6))
    fig_red, ax_red = plt.subplots(1, 1, figsize=(8, 6))

    # Interpolate for each redshift
    for z in input_redshifts:
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

    ax_blue.set_title(f"Blue, mass cut = {M_limit_effective:.1e}\nMass bins: {num_mass_bins}, z bins: {num_z_bins}")
    ax_red.set_title(f"Red, mass cut = {M_limit_effective:.1e}\nMass bins: {num_mass_bins}, z bins: {num_z_bins}")

    for ax in {ax_blue, ax_red}:
        ax.set_xlabel(r"Halo/Subhalo Mass ($M_{\odot}/h$)")
        ax.set_ylabel("Galaxy's Absolute Magnitude")

        ax.invert_yaxis()

        ax.set_xscale("log")

        ax.set_ylim([-18, -24])

        ax.legend(loc="lower right")

        ax.grid()

    fig_blue.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/blue_interpolation_byredshift_{pinocchio_particle_count}_{run_id}.pdf")
    fig_red.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/red_interpolation_byredshift_{pinocchio_particle_count}_{run_id}.pdf")

    #plt.show()

    # Close figures to reduce memory usage,
    # since pyplot keeps figures in memory until
    # the end of the program by default
    plt.close(fig_blue)
    plt.close(fig_red)
