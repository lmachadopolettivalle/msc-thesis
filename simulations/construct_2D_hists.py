# Description of the program:
# Using the subhalos extracted from PINOCCHIO,
# create a 2D histogram binned in mass and redshift,
# which will be used during the SHAM process.
# At the end, saves 2 files with 2D histograms,
# one for halos and one for subhalos.

# Author: Pascale Berner
# Co-Author: Luis Machado
# first written: 25.11.2022
# last adapted: 25.11.2022

# ------------------------------------------
# IMPORTS
# ------------------------------------------


print("Importing required libraries...")

import numpy as np
import os
import re
from tqdm import tqdm

print("Done importing libraries.")

# ------------------------------
# Interpolation parameters
# Can be modified and fine-tuned
# ------------------------------
num_z_bins = 150
num_mass_bins = 30


# ------------------------------
# Configurations, filenames and directories
# ------------------------------

# Cube root of number of particles used.
# This is present in the paths to different input files used in this script.
particle_count_pinocchio = 2048

pinocchio_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_{particle_count_pinocchio}" # Path to SLURM output from PINOCCHIO, which contains many useful details on the run

# Directory where subhalo files are stored
# These were generated with the subhalo code
dirname = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/"
halo_subhalo_files = "pinocchio_masked_halos_subhalos_plc"

# Obtain all subhalo files from the given directory
FILENAMES = [
    name for name in os.listdir(dirname)
    if halo_subhalo_files in name
]

hist_2D_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/halo_subhalo_plc/"
hist_2D_halo_file = "TESTpinocchio_masked_halos_hist2D"
hist_2D_subhalo_file = "TESTpinocchio_masked_subhalos_hist2D"


# ---------------------------------
# NO NEED TO MODIFY BELOW THIS LINE
# ---------------------------------

# Read following information from PINOCCHIO output file,
# after running PINOCCHIO
def read_pinocchio_config_details():
    with open(pinocchio_output_filename, 'r') as f:
        text = f.read()

    m_part = float(
        re.search("Particle Mass \(Msun/h\)\s+(.+)\n", text).groups()[0]
    )
    min_num_part = int(
        re.search("MinHaloMass \(particles\)\s+(\d+)\n", text).groups()[0]
    )
    z_max, z_min = (
        float(i)
        for i in re.search("The Past Light Cone will be reconstruct from z=(.+) to z=(.+)\n", text).groups()
    )

    return (
        m_part,
        min_num_part,
        z_min,
        z_max,
    )

m_part, min_num_part, z_min, z_max = read_pinocchio_config_details()

# Create bins for 2D histogram
min_mass = m_part * min_num_part # Msun/h
max_mass = 6.0e15 # Msun/h, hopefully high enough

bin_edges_mass = np.logspace(np.log10(min_mass), np.log10(max_mass), num=(num_mass_bins+1))
bin_edges_z = np.linspace(z_min, z_max, (num_z_bins+1))

hist_z_mass_halos, bin_edges_z, bin_edges_mass = np.histogram2d([], [], bins=(bin_edges_z, bin_edges_mass))
hist_z_mass_subs, bin_edges_z, bin_edges_mass = np.histogram2d([], [], bins=(bin_edges_z, bin_edges_mass))

# For each subhalo file, generate the corresponding histograms
print("Starting to process halo/subhalo files to create 2D histograms...")
for filename in tqdm(FILENAMES):
    print(f"Processing file {filename}...")
    masses, redshifts, is_halo = np.loadtxt(f"{dirname}/{filename}", unpack=True, usecols=(1, 2, 8))

    halo_mask = (is_halo == 1)
    subhalo_mask = (is_halo == 0)

    halo_masses = masses[halo_mask]
    halo_redshifts = redshifts[halo_mask]

    subhalo_masses = masses[subhalo_mask]
    subhalo_redshifts = redshifts[subhalo_mask]

    # ------------------------------
    # CALCULATE 2D HISTOGRAMS
    # ------------------------------
    hist_z_mass_halos_temp, bin_edges_z, bin_edges_mass = np.histogram2d(halo_redshifts, halo_masses, bins=(bin_edges_z, bin_edges_mass))
    hist_z_mass_subs_temp, bin_edges_z, bin_edges_mass = np.histogram2d(subhalo_redshifts, subhalo_masses, bins=(bin_edges_z, bin_edges_mass))

    hist_z_mass_halos += hist_z_mass_halos_temp
    hist_z_mass_subs += hist_z_mass_subs_temp

# ------------------------------
# SAVE 2D HISTOGRAMS TO FILES
# ------------------------------
output_2Dhist_halo = hist_2D_dir + hist_2D_halo_file # Note that savez adds a '.npz' extension
output_2Dhist_subhalo = hist_2D_dir + hist_2D_subhalo_file # Note that savez adds a '.npz' extension

np.savez(output_2Dhist_halo, hist_z_mass_halos=hist_z_mass_halos, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
np.savez(output_2Dhist_subhalo, hist_z_mass_subs=hist_z_mass_subs, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
