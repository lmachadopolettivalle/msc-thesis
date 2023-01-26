# Drawing Galaxies from Luminosity Functions

# UCat only, just get redshifts and absolute magnitudes, and red vs. blue
# No apparent magnitudes, no matching to PINOCCHIO
# Pixarea calculated from a mask

# Author: Pascale Berner
# Co-Author: Luis Machado
# first written: 09.11.2022
# last adapted: 21.11.2022
# partially copied from: sample_from_lum_fct_only.py

# ----------------------------------------------------
# IMPORTS
# -----------------------------------------------------
print("Importing required libraries...")

import h5py
import healpy as hp
import numpy as np
import os
import PyCosmo
import re
import scipy.integrate
from sham_model_constants import *
from ucat import galaxy_sampling_util, io_util
from ucat.galaxy_population_models import galaxy_luminosity_function

print("Done importing libraries.")

# ----------------------------------------------------
# PARAMETERS TO BE SET
# -----------------------------------------------------

# Configuration and cosmology parameters
m_max = -12

# Cube root of number of particles used.
# This is present in the paths to different input files used in this script.
particle_count_pinocchio = 2048

# File path names
pinocchio_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_deep_{particle_count_pinocchio}" # Path to SLURM output from PINOCCHIO, which contains many useful details on the run

outfile_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed_z1.5/outputs_sampling/"

outfile_name = "sampled_BASS"

outfile_z = f"{outfile_name}_z.npy"
outfile_absmag = f"{outfile_name}_absmag.npy"
outfile_redblue = f"{outfile_name}_redblue.npy"

# Path to file containing mask array as a HEALPix map
# If no mask is desired, set the filename to None,
# in which case the value DEFAULT_PIXAREA will be used,
# without any further declination cuts
infile_footprint = "/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy"
infile_footprint = None

NEST = True

MASK_CUTOFF = 0

# Set a default pixarea, which will be used in case there is no valid input mask
DEFAULT_PIXAREA = 0.007596 * 4 * np.pi # 10% of the sky

# Additionally, support a range of declinations
# In this case, the mask is further filtered to only keep its pixels within the
# provided declination range
# If this feature is not desired, set dec_min to -90 and dec_max to 90
dec_min = 30
dec_max = 90

# ----------------------------------------------------
# NO NEED TO MODIFY BELOW THIS LINE.
# ----------------------------------------------------

 # Create necessary directories
if os.path.isdir(outfile_dir):
    print(f"{outfile_dir} directory already exists.")
else:
    print(f"Creating new output directory, {outfile_dir} ...")
    os.mkdir(outfile_dir)
    print("Created output directory successfully.")

# ----------------------------------------------------
# LOAD MASK, CALCULATE PIXAREA
# -----------------------------------------------------
# Load pixel mask, if any is provided
pixarea = DEFAULT_PIXAREA

if infile_footprint is not None:
    print(f"Loading mask from {infile_footprint}...")

    if infile_footprint.endswith(".fit"):
        m_footprint = hp.fitsfunc.read_map(infile_footprint, nest=NEST)
    elif infile_footprint.endswith(".npy"):
        with open(infile_footprint, "rb") as f:
            m_footprint = np.load(f)
    else:
        print(f"WARNING: Invalid file format for HEALPix mask file {infile_footprint}. Using default pixarea {DEFAULT_PIXAREA:.2f} instead.")
        m_footprint = None

    if m_footprint is not None:
        NSIDE = hp.get_nside(m_footprint)
        NPIX = hp.nside2npix(NSIDE)

        # Use declination range to filter out mask further

        # Convert cutoff from degrees to colatitude
        # NOTE we flip min and max to higher and lower
        # due to the definition of colatitude
        higher_colatitude = np.pi/2 - np.radians(dec_min)
        lower_colatitude = np.pi/2 - np.radians(dec_max)

        # North
        # NOTE that query_strip has a bug using nest=True.
        # So if NEST is set to True above,
        # we first obtain the RING pixels and then convert to NEST manually.
        ring_dec_range_pixel_ids = hp.query_strip(
            NSIDE,
            theta1=lower_colatitude,
            theta2=higher_colatitude,
            inclusive=False,
            nest=False,
        )

        dec_range_pixel_ids = ring_dec_range_pixel_ids

        if NEST is True:
            ring_dec_range_map = np.zeros(NPIX)
            ring_dec_range_map[ring_dec_range_pixel_ids] = 1

            nest_dec_range_map = hp.reorder(ring_dec_range_map, r2n=True)
            nest_dec_range_pixel_ids = np.where(nest_dec_range_map > 0)[0]

            dec_range_pixel_ids = nest_dec_range_pixel_ids

        # Determine pixarea from provided map
        individual_pixarea = hp.pixelfunc.nside2pixarea(NSIDE, degrees=False)
        mask_pixel_ids = np.where(m_footprint > MASK_CUTOFF)[0]
        filtered_pixel_ids = np.intersect1d(mask_pixel_ids, dec_range_pixel_ids)

        pixarea = individual_pixarea * len(filtered_pixel_ids)
        print(f"Finished loading mask. Using pixarea {pixarea:.2f}, i.e. {100 * pixarea / (4 * np.pi):.2f}% of the sky.")
else:
    print(f"No mask file provided. Will use default pixarea {DEFAULT_PIXAREA:.2f} instead.")

# ----------------------------------------------------
# SET PYCOSMO COSMOLOGY
# -----------------------------------------------------

# Read following information from PINOCCHIO output file,
# after running PINOCCHIO
def read_pinocchio_config_details():
    with open(pinocchio_output_filename, 'r') as f:
        text = f.read()

    z_max, _ = (
        float(i)
        for i in re.search("The Past Light Cone will be reconstruct from z=(.+) to z=(.+)\n", text).groups()
    )
    omega_m = float(
        re.search("Omega0\s+(.+)\n", text).groups()[0]
    )
    omega_l = float(
        re.search("OmegaLambda\s+(.+)\n", text).groups()[0]
    )
    omega_b = float(
        re.search("OmegaBaryon\s+(.+)\n", text).groups()[0]
    )
    n_scalar = float(
        re.search("PrimordialIndex\s+(.+)\n", text).groups()[0]
    )
    sigma8 = float(
        re.search("Sigma8\s+(.+)\n", text).groups()[0]
    )
    h = float(
        re.search("Hubble100\s+(.+)\n", text).groups()[0]
    )

    return (
        z_max,
        omega_m,
        omega_l,
        omega_b,
        n_scalar,
        sigma8,
        h,
    )

z_max, omega_m, omega_l, omega_b, n_scalar, sigma8, h = read_pinocchio_config_details()


cosmo = PyCosmo.Cosmo()
cosmo.set(
    h=h,
    omega_b=omega_b,
    omega_m=omega_m,
    omega_l_in=omega_l,
    pk_norm_type="sigma8",
    pk_norm=sigma8,
    n=n_scalar,
)
cosmo.print_params()

# ----------------------------------------------------
# SAMPLE GALAXIES
# -----------------------------------------------------
n_gal_cal_blue = galaxy_sampling_util.NumGalCalculator(z_max, m_max, lum_fct_alpha_blue, m_star_par_blue, phi_star_par_blue, cosmo, pixarea)
n_gal_cal_red = galaxy_sampling_util.NumGalCalculator(z_max, m_max, lum_fct_alpha_red, m_star_par_red, phi_star_par_red, cosmo, pixarea)

n_gal_blue = n_gal_cal_blue()
n_gal_red = n_gal_cal_red()

z_mabs_sampler_blue = galaxy_luminosity_function.RedshiftAbsMagSampler(lum_fct_z_res, z_max, lum_fct_m_res, m_max, lum_fct_alpha_blue, m_star_par_blue, phi_star_par_blue, cosmo)

z_mabs_sampler_red = galaxy_luminosity_function.RedshiftAbsMagSampler(lum_fct_z_res, z_max, lum_fct_m_res, m_max, lum_fct_alpha_red, m_star_par_red, phi_star_par_red, cosmo)

z_blue, abs_mag_blue = z_mabs_sampler_blue(n_gal_blue)
z_red, abs_mag_red = z_mabs_sampler_red(n_gal_red)

index_blue = np.full(len(z_blue), BLUE, dtype=int)
index_red = np.full(len(z_red), RED, dtype=int)

n_tot = n_gal_red + n_gal_blue
z_tot = np.concatenate((z_red, z_blue))
abs_mag_tot = np.concatenate((abs_mag_red, abs_mag_blue))
index_tot = np.concatenate((index_red, index_blue))

# ----------------------------------------------------
# SAVE GALAXIES TO NPY FILES
# -----------------------------------------------------
with open(outfile_dir + outfile_z, 'wb') as f:
    np.save(f, z_tot)
with open(outfile_dir + outfile_absmag, 'wb') as f:
    np.save(f, abs_mag_tot)
with open(outfile_dir + outfile_redblue, 'wb') as f:
    np.save(f, index_tot)

# ----------------------------------------------------
# SOME PRINT STATEMENTS
# -----------------------------------------------------

print('z_max = ' + str(z_max))
print('pixarea = ' + str(pixarea))
print('outfile_dir = ' + outfile_dir)
print('outfile_z = ' + outfile_z)
print('outfile_absmag = ' + outfile_absmag)
print('outfile_redblue = ' + outfile_redblue)
print('output: abs_mag_tot, z_tot, index_tot')

