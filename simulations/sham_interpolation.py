# Description of the program:
# read in sampled galaxies
# read in 2D histograms of halos / subhalos
# define 2D interpolation
# read in halo-subhalo catalog
# populate halos and subhalos with galaxies
# calculate apparent magnitudes, apply cuts, save output
# Note: with a mask, using multiple halo-subhalo plc files

# Author: Pascale Berner
# Co-Author: Luis Machado
# first written: 17.11.2022
# last adapted: 28.11.2022
# partially copied from: sample_from_lum_fct_interp_pin_desi.py etc.

# ----------------------------------------------------
# IMPORTS
# -----------------------------------------------------
print("Importing required libraries...")

import argparse
import concurrent.futures
import h5py
from manage_parameter_space import get_details_of_run
import numpy as np
import os
import pandas as pd
import PyCosmo
import re
from scipy.interpolate import griddata
import scipy.integrate
from sham_model_constants import *
from ucat import galaxy_sampling_util, io_util
from ucat.galaxy_population_models import galaxy_luminosity_function, galaxy_sed

print("Done importing libraries.")

# ----------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--run_id", type=int, required=True)
parser.add_argument("--particle_count_pinocchio", type=int, required=True)
parser.add_argument("--region", type=str, required=True)

args = parser.parse_args()

run_id = args.run_id
particle_count_pinocchio = args.particle_count_pinocchio
region = args.region # BASS or DECaLS

# Get details of run
run_details = get_details_of_run(run_id)
# SHAM parameters
M_limit = run_details["mass_cut"] # mass limit for assigning blue or red galaxies to halos, [Msun/h]
M_limit_effective = M_limit
quenching_time = run_details["quenching_time"]

# Desired filters
# NOTE that the filter lum_fct_filter_band (defined in the other .py script)
# is added to this list, to work well with the luminosity function
if region == "BASS":
    print("Using BASS/MzLS filters")
    desired_filters = {
        "g": "BASSMzLS_g",
        "r": "BASSMzLS_r",
        "z": "BASSMzLS_z",
    }
else:
    print("Using DECam filters")
    desired_filters = {
        "g": "DECam_g",
        "r": "DECam_r",
        "z": "DECam_z",
    }


pinocchio_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_{particle_count_pinocchio}" # Path to SLURM output from PINOCCHIO, which contains many useful details on the run

# Create output directory for SHAM results
outfile_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/{run_id}/interpolation_outputs/"
outfile_galaxies_base = "ucat_sorted_app_mag_interp_"
if os.path.isdir(outfile_dir):
    print(f"{outfile_dir} directory already exists.")
else:
    print(f"Creating new output directory, {outfile_dir} ...")
    os.mkdir(outfile_dir)
    print("Created output directory successfully.")

# Paths to hist2D files generated by subhalo code
infile_hist2D_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/{run_id}/2D_histograms/"
infile_hist_red = 'pinocchio_masked_red_hist2D.npz'
infile_hist_blue = 'pinocchio_masked_blue_hist2D.npz'

# Paths to halos and subhalos catalogs generated by PINOCCHIO + the subhalo code
infile_halos_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/"
infile_halos = 'pinocchio_masked_halos_subhalos_plc'

# Paths to UCat output files, generated by UCat sampling
# Directory where sampling outputs are stored
infile_ucat_sampling_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/outputs_sampling/"

infile_ucat_z = "sampled_BASS_z.npy"
infile_ucat_absmag = "sampled_BASS_absmag.npy"
infile_ucat_redblue = "sampled_BASS_redblue.npy"

# Output file names
output_interp_red = 'red_lim_interp.npz'
output_interp_blue = 'blue_lim_interp.npz'


# apparent magnitude limits, in the main band (typically i or r)
BAND_USED_FOR_MAG_CUTS = "r"
gals_mag_max = 19.5
gals_mag_min = 10


# ----------------------------------------------------
# NO NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------
# Add the luminosity function-related filter to the list of desired filters
loop_filter_names = [lum_fct_filter_band] + list(desired_filters.values())

# ----------------------------------------------------
# SET PYCOSMO COSMOLOGY
# -----------------------------------------------------
# Read following information from PINOCCHIO output file,
# after running PINOCCHIO
def read_pinocchio_config_details():
    with open(pinocchio_output_filename, 'r') as f:
        text = f.read()

    # num_files corresponds to NumFiles in PINOCCHIO output file
    # Note that it may differ from the requested NumFiles in the input parameter file. PINOCCHIO requires that NumFiles divides NumTasks, and may change the NumFiles value if not.
    num_files = int(
        re.search("NumFiles\s+(\d+)\n", text).groups()[0]
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

    return (
        num_files,
        omega_m,
        omega_l,
        omega_b,
        n_scalar,
        sigma8,
    )

num_files, omega_m, omega_l, omega_b, n_scalar, sigma8 = read_pinocchio_config_details()

if num_files == 1:
    file_ending = ['']
else:
    file_ending = [f'.{i}' for i in range(num_files)]


cosmo = PyCosmo.Cosmo()
cosmo.set(
    omega_b=omega_b,
    omega_m=omega_m,
    omega_l_in=omega_l,
    pk_norm_type="sigma8",
    pk_norm=sigma8,
    n=n_scalar,
)
cosmo.print_params()

# ----------------------------------------------------
# DEFINE UFIG CLASSES
# -----------------------------------------------------
class Filter(object):
	"""
	Class representing a filter defined on a grid of wavelengths.
	"""

	def __init__(self, lam, amp):
		"""
		Class initializer.
		:param lam: Wavelength grid the filter is defined on.
		:param amp: Amplitude of the filter at each grid point.
		"""
		self.lam = lam
		self.amp = amp

# MAGNITUDE CALCULATOR TABLE
class MagCalculatorTable(object):
	def __init__(self, filters):
		self.c = 3 * 1e8 * 1e6  # speed of light, units: micrometre/s
		self.z_grid, self.excess_b_v_grid = io_util.load_from_hdf5(templates_int_tables_file_name, ('z', 'E(B-V)'), root_path=maps_remote_dir)
		self.templates_int_tables_file_name = templates_int_tables_file_name
		self.maps_remote_dir = maps_remote_dir
		self.update_filters(filters)

	def __call__(self, redshifts, excess_b_v, coeffs, filter_names):
		magnitudes = {}
		z_ind = galaxy_luminosity_function.find_closest_ind(self.z_grid, redshifts)
		excess_b_v_ind = galaxy_luminosity_function.find_closest_ind(self.excess_b_v_grid, excess_b_v)

		for filter_name in filter_names:
			templates_int_tables = self.templates_int_table_dict[filter_name]
			mags = np.zeros(redshifts.size, dtype=np.float64)
			for i in range(coeffs.shape[1]):
				mags += coeffs[:, i] * templates_int_tables[i][z_ind, excess_b_v_ind]
			mags = -2.5 * np.log10(mags) - 48.6
			mags[np.isnan(mags)] = np.inf
			magnitudes[filter_name] = mags

		return magnitudes

	def update_filters(self, filters):
		self.get_n_templates(filters.amp.keys())
		self.templates_int_table_dict = {}
		self.filter_norm_dict = {}

		for filter_name in filters.amp.keys():
			self.templates_int_table_dict[filter_name] = \
				io_util.load_from_hdf5(self.templates_int_tables_file_name,
				['integrals/{}/template_{}'.format(filter_name, i) for i in
				range(self.n_templates_dict[filter_name])], root_path=self.maps_remote_dir)
			self.filter_norm_dict[filter_name] = scipy.integrate.simps(filters.amp[filter_name] / filters.lam[filter_name],
																		x=filters.lam[filter_name])

	def get_n_templates(self, filter_names):
		with h5py.File(io_util.get_abs_path(self.templates_int_tables_file_name, root_path=self.maps_remote_dir), mode='r') as f:
			self.n_templates_dict = {filter_name: len(list(filter(lambda k: k.startswith('template_'),
									f['integrals/{}'.format(filter_name)].keys())))
									for filter_name in filter_names}

# ----------------------------------------------------–
# DEFINE FUNCTIONS
# -----------------------------------------------------
# DEFINE MAGNITUDE CALCULATOR
print("Starting magnitude calculations...")

MAGNITUDES_CALCULATOR = {'table': MagCalculatorTable}
filter_wavelengths = io_util.load_from_hdf5(filters_file_name, [i+'/lam' for i in loop_filter_names], root_path=maps_remote_dir)
print("Loaded first hdf5")
filter_amplitudes = io_util.load_from_hdf5(filters_file_name, [i+'/amp' for i in loop_filter_names], root_path=maps_remote_dir)
print("Loaded second hdf5")
filter_lam = {loop_filter_names[i]: filter_wavelengths[i] for i in range(len(loop_filter_names))}
filter_amp = {loop_filter_names[i]: filter_amplitudes[i] for i in range(len(loop_filter_names))}

loop_filters = Filter(lam=filter_lam, amp=filter_amp)
mag_calc = MAGNITUDES_CALCULATOR[magnitude_calculation](loop_filters)
n_templates = mag_calc.n_templates_dict[lum_fct_filter_band]
print("Finished n_templates")

# ----------------------------------------------------–
# LOAD GALAXIES
# -----------------------------------------------------

abs_mag = np.load(infile_ucat_sampling_dir + infile_ucat_absmag)
z_ucat = np.load(infile_ucat_sampling_dir + infile_ucat_z)
blue_red = np.load(infile_ucat_sampling_dir + infile_ucat_redblue)

# Sort by absolute magnitude before performing SHAM
abs_mag_inds_sorted = np.argsort(abs_mag)
abs_mag = abs_mag[abs_mag_inds_sorted]
z_ucat = z_ucat[abs_mag_inds_sorted]
blue_red = blue_red[abs_mag_inds_sorted]

# Split between red and blue galaxies
abs_mag_red = abs_mag[blue_red == RED]
z_red = z_ucat[blue_red == RED]
abs_mag_blue = abs_mag[blue_red == BLUE]
z_blue = z_ucat[blue_red == BLUE]

# ----------------------------------------------------–
# LOAD HALO-SUBHALO HISTOGRAMS
# -----------------------------------------------------
print("Loading hist reds")
with np.load(infile_hist2D_dir + infile_hist_red) as data:
    hist_z_mass_red=data['hist_z_mass_red']
    bin_edges_z=data['bin_edges_z']
    bin_edges_mass=data['bin_edges_mass']

print("Loading hist blues")
with np.load(infile_hist2D_dir + infile_hist_blue) as data:
    hist_z_mass_blue=data['hist_z_mass_blue']

num_z_bins = len(bin_edges_z) - 1
num_mass_bins = len(bin_edges_mass) - 1

# ----------------------------------------------------–
# CREATE 2D HISTOGRAMS FOR RED vs. BLUE
# -----------------------------------------------------
# NOTE: This part changes, in case the SHAM method changes.

# ----------------------------------------------------–
# CALCULATE 2D ARRAYS FOR INTERPOLATION
# -----------------------------------------------------
# empty arrays for values; limit meaning at the edges of 2D histogram
lim_abs_mag_red = np.zeros((num_z_bins+1, num_mass_bins+1))
lim_abs_mag_blue = np.zeros((num_z_bins+1, num_mass_bins+1))
lim_abs_mag_red[:] = np.nan
lim_abs_mag_blue[:] = np.nan

num_z = num_z_bins + 1
num_mass = num_mass_bins + 1

print("Starting loop redshift")
for i in range(num_z - 1):  # loop over redshift, starting at low z
	abs_mag_red_i = abs_mag_red[(z_red > bin_edges_z[i]) & (z_red <= bin_edges_z[i+1])]
	abs_mag_blue_i = abs_mag_blue[(z_blue > bin_edges_z[i]) & (z_blue <= bin_edges_z[i+1])]

	ind_red = 0
	ind_blue = 0

	for j in reversed(range(num_mass - 1)):  # loop over mass bins, starting at hight mass
		num_blue = hist_z_mass_blue[i,j]
		if int(ind_blue + num_blue + 1) >= len(abs_mag_blue_i):  # ensure there are enough galaxies
			break
		if num_blue > 0:
			#mean_abs_mag_blue[i,j] = np.mean(abs_mag_blue_i[ind_blue : int(ind_blue + num_blue)])
			lim_abs_mag_blue[i+1,j] = abs_mag_blue_i[int(ind_blue + num_blue)]
			if np.isnan(lim_abs_mag_blue[i,j]):  # check if the neighboring bin in redshift is empty
				lim_abs_mag_blue[i,j] = abs_mag_blue_i[int(ind_blue + num_blue)]
			if np.isnan(lim_abs_mag_blue[i+1,j+1]):  # check if the neighboring bin in mass is empty
				lim_abs_mag_blue[i+1,j+1] = abs_mag_blue_i[ind_blue]
			ind_blue =  int(ind_blue + num_blue)

	for j in reversed(range(num_mass - 1)):  # loop over mass bins, starting at hight mass
		num_red = hist_z_mass_red[i,j]
		if int(ind_red + num_red + 1) >= len(abs_mag_red_i):  # ensure there are enough galaxies
			break
		if num_red > 0:
			#mean_abs_mag_red[i,j] = np.mean(abs_mag_red_i[ind_red : int(ind_red + num_red)])
			lim_abs_mag_red[i+1,j] = abs_mag_red_i[int(ind_red + num_red)]
			if np.isnan(lim_abs_mag_red[i,j]):  # check if the neighboring bin in redshift is empty
				lim_abs_mag_red[i,j] = abs_mag_red_i[int(ind_red + num_red)]
			if np.isnan(lim_abs_mag_red[i+1,j+1]):  # check if the neighboring bin in mass is empty
				lim_abs_mag_red[i+1,j+1] = abs_mag_red_i[ind_red]
			ind_red = int(ind_red + num_red)

# ----------------------------------------------------–
# SAVE INTERPOLATION PROPERTIES
# -----------------------------------------------------
print("Saving npz files")
np.savez(outfile_dir + output_interp_red, lim_abs_mag_red=lim_abs_mag_red, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
np.savez(outfile_dir + output_interp_blue, lim_abs_mag_blue=lim_abs_mag_blue, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)

# ----------------------------------------------------–
# PREPARE STACKED MASS EDGES
# -----------------------------------------------------
mass_edges_stacked = np.reshape(np.tile(bin_edges_mass, len(bin_edges_z)), (len(bin_edges_mass)*len(bin_edges_z),))
z_edges_stacked = np.reshape(np.ndarray.flatten(np.transpose(np.tile(bin_edges_z, (len(bin_edges_mass), 1)))), (len(bin_edges_mass)*len(bin_edges_z),))

# ----------------------------------------------------–
# LOAD HALO-SUBHALO FILES INTO MEMORY
# -----------------------------------------------------

def process_halo_subhalo_file(i):
    """Given a file index i (int), process the corresponding halo subhalo file,
    and generate a catalog of galaxies assigned to the loaded halos and subhalos
    using SHAM.
    """
    # ----------------------------------------------------–
    # LOAD HALO-SUBHALO FILE
    # -----------------------------------------------------
    filename = infile_halos_dir + infile_halos + file_ending[i] + '.txt'
    data = pd.read_csv(filename, sep='\s+', lineterminator='\n', header=None, index_col=None, skipinitialspace=True).values

    mass = data[:, 1].copy()
    z_pin = data[:, 2].copy()
    x_coord_pin = data[:, 3].copy()
    y_coord_pin = data[:, 4].copy()
    z_coord_pin = data[:, 5].copy()
    host_sub = data[:, 8].copy()
    delta_t = data[:, 9].copy()

    # ----------------------------------------------------–
    # DIVIDE HALOS AND SUBHALOS INTO RED AND BLUE
    # -----------------------------------------------------
    mask_red = ((host_sub == 1) & (mass > M_limit)) | ((host_sub == 0) & (delta_t > quenching_time))
    mask_blue = ~ mask_red

    n_blue = len(z_pin[mask_blue])

    n_uncut_temp = len(z_pin)
    n_blue_uncut_temp = n_blue
    # ----------------------------------------------------–
    # GET ABSOLUTE MAGNITUDES FOR NEW GALAXIES
    # -----------------------------------------------------
    abs_mag_red_pin = griddata((z_edges_stacked, np.log10(mass_edges_stacked)), np.ndarray.flatten(lim_abs_mag_red), (z_pin[mask_red], np.log10(mass[mask_red])))
    abs_mag_blue_pin = griddata((z_edges_stacked, np.log10(mass_edges_stacked)), np.ndarray.flatten(lim_abs_mag_blue), (z_pin[mask_blue], np.log10(mass[mask_blue])))
    # ----------------------------------------------------–
    # APPEND RED AND BLUE GALAXIES
    # -----------------------------------------------------
    z_temp = np.append(z_pin[mask_blue], z_pin[mask_red])
    abs_mag_temp = np.append(abs_mag_blue_pin, abs_mag_red_pin)
    blue_red_temp = np.ones_like(z_temp, dtype=np.int16)
    blue_red_temp = np.full(len(z_temp), BLUE, dtype=int)
    blue_red_temp[n_blue:] = RED
    halo_mass_temp = np.append(mass[mask_blue], mass[mask_red])
    x_coord_temp = np.append(x_coord_pin[mask_blue], x_coord_pin[mask_red])
    y_coord_temp = np.append(y_coord_pin[mask_blue], y_coord_pin[mask_red])
    z_coord_temp = np.append(z_coord_pin[mask_blue], z_coord_pin[mask_red])
    host_sub_index_temp = np.append(host_sub[mask_blue], host_sub[mask_red])
    time_since_merger_temp = np.append(delta_t[mask_blue], delta_t[mask_red])
    # ----------------------------------------------------–
    # CALCULATE TEMPLATE COEFFICIENTS
    # -----------------------------------------------------
    # Draw template coefficients
    template_coeffs = np.empty((len(z_temp), n_templates))
    template_coeffs[:n_blue] = galaxy_sed.sample_template_coeff_dirichlet(z_pin[mask_blue], template_coeff_alpha0_blue, template_coeff_alpha1_blue, template_coeff_z1_blue, template_coeff_weight_blue)
    template_coeffs[n_blue:] = galaxy_sed.sample_template_coeff_dirichlet(z_pin[mask_red], template_coeff_alpha0_red, template_coeff_alpha1_red, template_coeff_z1_red, template_coeff_weight_red)
    # Calculate absolute magnitudes according to coefficients and adjust coefficients according to drawn magnitudes
    template_coeffs *= np.expand_dims(10 ** (0.4 * (mag_calc(np.zeros_like(z_temp), np.zeros_like(z_temp), template_coeffs, [lum_fct_filter_band])[lum_fct_filter_band] - abs_mag_temp)), -1)
    # Transform to apparent coefficients
    lum_dist = galaxy_sampling_util.apply_pycosmo_distfun(cosmo.background.dist_lum_a, z_temp)
    template_coeffs *= np.expand_dims((10e-6 / lum_dist) ** 2 / (1 + z_temp), -1)
    # ----------------------------------------------------–
    # SET EXTINCTION MAP
    # -----------------------------------------------------
    # I'm setting them to zero, since I'm ignoring extinction values, as I don't have pixel coordinates
    excess_b_v = np.zeros(len(z_temp))
    # ----------------------------------------------------–
    # CALCULATE APPARENT MAGNITUDES
    # -----------------------------------------------------
    temp_app_mag_dict = {}
    for k in app_mag_dict.keys():
        temp_app_mag_dict[k] = mag_calc(
            z_temp,
            excess_b_v,
            template_coeffs,
            [desired_filters[k]]
        )[desired_filters[k]]

    # ----------------------------------------------------–
    # APPLY CUTS IN APPARENT MAGNITUDES
    # -----------------------------------------------------
    mask_mag_range = (temp_app_mag_dict[BAND_USED_FOR_MAG_CUTS] >= gals_mag_min) & (temp_app_mag_dict[BAND_USED_FOR_MAG_CUTS] <= gals_mag_max)

    z_temp = z_temp[mask_mag_range]
    abs_mag_temp = abs_mag_temp[mask_mag_range]
    blue_red_temp = blue_red_temp[mask_mag_range]

    for k in temp_app_mag_dict.keys():
        temp_app_mag_dict[k] = temp_app_mag_dict[k][mask_mag_range]

    halo_mass_temp = halo_mass_temp[mask_mag_range]
    x_coord_temp = x_coord_temp[mask_mag_range]
    y_coord_temp = y_coord_temp[mask_mag_range]
    z_coord_temp = z_coord_temp[mask_mag_range]
    host_sub_index_temp = host_sub_index_temp[mask_mag_range]
    time_since_merger_temp = time_since_merger_temp[mask_mag_range]

    return (
        n_uncut_temp,
        n_blue_uncut_temp,
        z_temp,
        abs_mag_temp,
        blue_red_temp,
        temp_app_mag_dict,
        halo_mass_temp,
        x_coord_temp,
        y_coord_temp,
        z_coord_temp,
        host_sub_index_temp,
        time_since_merger_temp,
    )

# ----------------------------------------------------–
# PREPARE EMPTY ARRAYS FOR GALAXIES
# -----------------------------------------------------
n_uncut = 0
n_blue_uncut = 0

z = np.array([])
abs_mag = np.array([])
blue_red = np.array([])

app_mag_dict = {
    k: np.array([])
    for k in desired_filters.keys()
}

halo_mass = np.array([])
x_coord = np.array([])
y_coord = np.array([])
z_coord = np.array([])
host_sub_index = np.array([])
time_since_merger = np.array([])

print("Processing halo subhalo files...")
for i in range(num_files):
    print(f"Processing file index {i}...")

    n_uncut_temp, n_blue_uncut_temp, z_temp, abs_mag_temp, blue_red_temp, temp_app_mag_dict, halo_mass_temp, x_coord_temp, y_coord_temp, z_coord_temp, host_sub_index_temp, time_since_merger_temp = process_halo_subhalo_file(i)

    n_uncut += n_uncut_temp
    n_blue_uncut += n_blue_uncut_temp

    z = np.append(z, z_temp)
    abs_mag = np.append(abs_mag, abs_mag_temp)
    blue_red = np.append(blue_red, blue_red_temp)

    for k in app_mag_dict.keys():
        app_mag_dict[k] = np.append(app_mag_dict[k], temp_app_mag_dict[k])

    halo_mass = np.append(halo_mass, halo_mass_temp)
    x_coord = np.append(x_coord, x_coord_temp)
    y_coord = np.append(y_coord, y_coord_temp)
    z_coord = np.append(z_coord, z_coord_temp)
    host_sub_index = np.append(host_sub_index, host_sub_index_temp)
    time_since_merger = np.append(time_since_merger, time_since_merger_temp)

# ----------------------------------------------------–
# SAVE OUTPUT
# -----------------------------------------------------
# each property in a separate npy file
with open(outfile_dir + outfile_galaxies_base + "z.npy", 'wb') as f:
    np.save(f, z)
with open(outfile_dir + outfile_galaxies_base + "abs_mag.npy", 'wb') as f:
    np.save(f, abs_mag)
with open(outfile_dir + outfile_galaxies_base + "blue_red.npy", 'wb') as f:
    np.save(f, blue_red)

for k, values in app_mag_dict.items():
    with open(outfile_dir + outfile_galaxies_base + f"app_mag_{k}.npy", 'wb') as f:
        np.save(f, values)

with open(outfile_dir + outfile_galaxies_base + "halo_mass.npy", 'wb') as f:
    np.save(f, halo_mass)
with open(outfile_dir + outfile_galaxies_base + "x_coord.npy", 'wb') as f:
    np.save(f, x_coord)
with open(outfile_dir + outfile_galaxies_base + "y_coord.npy", 'wb') as f:
    np.save(f, y_coord)
with open(outfile_dir + outfile_galaxies_base + "z_coord.npy", 'wb') as f:
    np.save(f, z_coord)
with open(outfile_dir + outfile_galaxies_base + "host_sub_index.npy", 'wb') as f:
    np.save(f, host_sub_index)
with open(outfile_dir + outfile_galaxies_base + "time_since_merger.npy", 'wb') as f:
    np.save(f, time_since_merger)

# ----------------------------------------------------–
# SOME PRINT STATEMENTS
# -----------------------------------------------------

print('Code run.')
print('num_files = ' + str(num_files))
print('infile_halos_dir = ' + infile_halos_dir)
print('infile_hist_red = ' + infile_hist_red)
print('infile_hist_blue = ' + infile_hist_blue)
print('infile_halos = ' + infile_halos)
print('infile_ucat_z = ' + infile_ucat_z)
print('infile_ucat_absmag = ' + infile_ucat_absmag)
print('infile_ucat_redblue = ' + infile_ucat_redblue)
print('output_interp_red = ' + output_interp_red)
print('output_interp_blue = ' + output_interp_blue)
print('outfile_galaxies_base = ' + outfile_galaxies_base)
print('M_limit = ' + str(M_limit))
print('M_limit_effective = ' + str(M_limit_effective))
print('gals_mag_max = ' + str(gals_mag_max))
print('gals_mag_min = ' + str(gals_mag_min))
print('omega_b = ' + str(omega_b))
print('omega_m = ' + str(omega_m))
print('omega_l_in = ' + str(omega_l))
print('pk_norm (sigma8) = ' + str(sigma8))
print('n = ' + str(n_scalar))
print('filters = ' + str([loop_filter_names[i] for i in range(len(loop_filter_names))]))
print('total number of galaxies without magnitude cuts = ' + str(n_uncut))
print('number of blue galaxies without magnitude cuts = ' + str(n_blue_uncut))
print('number of red galaxies without magnitude cuts = ' + str(n_uncut - n_blue_uncut))
print('total number of galaxies after magnitude cuts = ' + str(len(z)))
print('number of blue galaxies after magnitude cuts = ' + str(np.sum(blue_red)))
print('number of red galaxies after magnitude cuts = ' + str(len(z) - np.sum(blue_red)))
print('number of galaxies in halos after magnitude cuts = ' + str(np.sum(host_sub_index)))
print('number of galaxies in subhalos after magnitude cuts = ' + str(len(z) - np.sum(host_sub_index)))
