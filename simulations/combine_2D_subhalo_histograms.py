# During the subhalo construction code,
# we compute the halo and subhalo 2D histograms
# (binned in redshift and mass), and
# save them to different files.
# This script adds up all of the created 2D histograms,
# and saves the complete 2D histogram to a new file

import numpy as np
import os

# Cube root of number of particles used.
# This is present in the paths to different input files used in this script.
particle_count_pinocchio = 2048

hist_2D_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/halo_subhalo_plc/"

hist_2D_halo_file = "pinocchio_masked_halos_hist2D"
hist_2D_subhalo_file = "pinocchio_masked_subhalos_hist2D"

FILENAMES = sorted(os.listdir(hist_2D_dir))

# Loop through halo and subhalo files,
# accumulating the histograms
hist_z_mass = {
    "halos": None,
    "subs": None,
}

for filename in FILENAMES:
    print(f"Processing file {filename}...")
    # Should not include previously completed accumulation
    if (filename == f"{hist_2D_halo_file}.npz") or (filename == f"{hist_2D_subhalo_file}.npz"):
        print("Skipping file since it contains old, accumulated histogram data...")
        continue

    if hist_2D_halo_file in filename:
        halo_or_subhalo = "halos"
    elif hist_2D_subhalo_file in filename:
        halo_or_subhalo = "subs"
    else:
        print(f"{filename} does not match file format for either halos or subhalos.")
        raise ValueError

    with np.load(f"{hist_2D_dir}/{filename}") as data:
        hist_z_mass_temp = data[f"hist_z_mass_{halo_or_subhalo}"]
        bin_edges_z = data["bin_edges_z"]
        bin_edges_mass = data["bin_edges_mass"]

    if hist_z_mass[halo_or_subhalo] is None:
        hist_z_mass[halo_or_subhalo] = hist_z_mass_temp
    else:
        hist_z_mass[halo_or_subhalo] += hist_z_mass_temp

    # Remove temporary file after it has been added to the cumulative histogram
    #os.remove(f"{hist_2D_dir}/{filename}")

# -----------------------------------------------------
# SAVE 2D HISTOGRAMS
# -----------------------------------------------------
if( hist_z_mass["halos"] is None) or (hist_z_mass["subs"] is None):
    print("No histograms have been processed. Will not write output data.")
else:
    np.savez(hist_2D_dir + hist_2D_halo_file, hist_z_mass_halos=hist_z_mass["halos"], bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
    np.savez(hist_2D_dir + hist_2D_subhalo_file, hist_z_mass_subs=hist_z_mass["subs"], bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
