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

import argparse
import concurrent.futures
import numpy as np
import os
import pandas as pd
import re

from manage_parameter_space import get_details_of_run

print("Done importing libraries.")

# ------------------------------
# Interpolation parameters
# Can be modified and fine-tuned
# ------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--run_id", type=int, required=True)
parser.add_argument("--particle_count_pinocchio", type=int, required=True)

args = parser.parse_args()

run_id = args.run_id
particle_count_pinocchio = args.particle_count_pinocchio

# Get details of run
run_details = get_details_of_run(run_id)
num_z_bins = run_details["num_z_bins"]
num_mass_bins = run_details["num_mass_bins"]
M_limit = run_details["mass_cut"]
quenching_time = run_details["quenching_time"]

# ------------------------------
# Configurations, filenames and directories
# ------------------------------

pinocchio_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_deep_{particle_count_pinocchio}" # Path to SLURM output from PINOCCHIO, which contains many useful details on the run

# Directory where subhalo files are stored
# These were generated with the subhalo code
dirname = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed_z1.5/"
halo_subhalo_files = "pinocchio_masked_halos_subhalos_plc"

# Obtain all subhalo files from the given directory
FILENAMES = [
    name for name in os.listdir(dirname)
    if halo_subhalo_files in name
]

run_directory = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed_z1.5/{run_id}/"
if os.path.isdir(run_directory):
    print(f"{run_directory} directory already exists.")
else:
    print(f"Creating new output directory, {run_directory} ...")
    os.mkdir(run_directory)
    print("Created output directory successfully.")

hist_2D_dir = f"{run_directory}/2D_histograms/"

hist_2D_red_file = "pinocchio_masked_red_hist2D"
hist_2D_blue_file = "pinocchio_masked_blue_hist2D"

if os.path.isdir(hist_2D_dir):
    print(f"{hist_2D_dir} directory already exists.")
else:
    print(f"Creating new output directory, {hist_2D_dir} ...")
    os.mkdir(hist_2D_dir)
    print("Created output directory successfully.")

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

hist_z_mass_red, bin_edges_z, bin_edges_mass = np.histogram2d([], [], bins=(bin_edges_z, bin_edges_mass))
hist_z_mass_blue, bin_edges_z, bin_edges_mass = np.histogram2d([], [], bins=(bin_edges_z, bin_edges_mass))

# For each subhalo file, generate the corresponding histograms
print("Starting to process halo/subhalo files to create 2D histograms...")
def process_halo_subhalo_file(filename):
    global bin_edges_z
    global bin_edges_mass
    print(f"Processing file {filename}...")
    data = pd.read_csv(f"{dirname}/{filename}", sep='\s+', lineterminator='\n', header=None, index_col=None, skipinitialspace=True).values
    masses = data[:, 1]
    redshifts = data[:, 2]
    is_halo = data[:, 8]
    time_since_merger = data[:, 9]

    red_mask = ((is_halo == 1) & (masses > M_limit)) | ((is_halo == 0) & (time_since_merger > quenching_time))
    blue_mask = ((is_halo == 1) & (masses <= M_limit)) | ((is_halo == 0) & (time_since_merger <= quenching_time))

    red_masses = masses[red_mask]
    red_redshifts = redshifts[red_mask]

    blue_masses = masses[blue_mask]
    blue_redshifts = redshifts[blue_mask]

    # ------------------------------
    # CALCULATE 2D HISTOGRAMS
    # ------------------------------
    hist_z_mass_red_temp, bin_edges_z_temp, bin_edges_mass_temp = np.histogram2d(red_redshifts, red_masses, bins=(bin_edges_z, bin_edges_mass))
    hist_z_mass_blue_temp, bin_edges_z_temp, bin_edges_mass_temp = np.histogram2d(blue_redshifts, blue_masses, bins=(bin_edges_z, bin_edges_mass))

    return (hist_z_mass_red_temp, hist_z_mass_blue_temp, bin_edges_z_temp, bin_edges_mass_temp)


with concurrent.futures.ProcessPoolExecutor() as executor:
    for result in executor.map(process_halo_subhalo_file, FILENAMES):
        hist_z_mass_red_temp, hist_z_mass_blue_temp, bin_edges_z, bin_edges_mass = result

        hist_z_mass_red += hist_z_mass_red_temp
        hist_z_mass_blue += hist_z_mass_blue_temp

# ------------------------------
# SAVE 2D HISTOGRAMS TO FILES
# ------------------------------
output_2Dhist_red = hist_2D_dir + hist_2D_red_file # Note that savez adds a '.npz' extension
output_2Dhist_blue = hist_2D_dir + hist_2D_blue_file # Note that savez adds a '.npz' extension

np.savez(output_2Dhist_red, hist_z_mass_red=hist_z_mass_red, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
np.savez(output_2Dhist_blue, hist_z_mass_blue=hist_z_mass_blue, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
