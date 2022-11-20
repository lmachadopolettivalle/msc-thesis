# Description of the program:
# read in sampled galaxies
# read in 2D histograms of halos / subhalos
# define 2D interpolation
# read in halo-subhalo catalog
# populate halos and subhalos with galaxies
# calculate apparent magnitudes, apply cuts, save output
# Note: with a mask, using multiple halo-subhalo plc files

# Author: Pascale Berner
# first written: 17.11.2022
# last adapted: 18.11.2022
# partially copied from: sample_from_lum_fct_interp_pin_desi.py etc.

# ----------------------------------------------------
# IMPOTRS
# -----------------------------------------------------
from ucat import galaxy_sampling_util
from ucat import io_util
from ucat.galaxy_population_models import galaxy_luminosity_function
from ucat.galaxy_population_models import galaxy_sed
from scipy.interpolate import griddata
import h5py
import scipy.integrate
import PyCosmo
import numpy as np

# ----------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------
dirname = './' # recommendation: create a directory for each simulation, since many files are saved
num_files = 24
# only one file: file_ending = ['']
infile_halos_dir = '/cluster/scratch/bernerp/pinocchio/pin_500mpc_2048_11/halo_subhalo_plc/'
infile_hist_halos = 'pinocchio_masked_halos_hist2D'
infile_hist_subhalos = 'pinocchio_masked_subhalos_hist2D'
infile_halos = 'pinocchio_masked_halos_subhalos_plc'

infile_ucat_z = 'sampled_z_desy1_20221118_1.npy'
infile_ucat_absmag = 'sampled_absmag_desy1_20221118_1.npy'
infile_ucat_redblue = 'sampled_redblue_desy1_20221118_1.npy'

output_interp_red = 'lim_red_interp_20221118_1.npz'
output_interp_blue = 'lim_red_interp_20221118_1.npz'

outfile_galaxies_base = '/cluster/scratch/bernerp/ucat/des/ucat_sorted_app_mag_interp_20221118_1_'

# SHAM parameters
M_limit = 8.e12 # mass limit for assigning blue or red galaxies to halos, [Msun/h]

# apparent magnitude limits, in the main band (typically i or r)
gals_mag_max = 20.
gals_mag_min = 10.

if num_files == 1:
    file_ending = ['']
else:
    file_ending = [f'.{i}' for i in range(num_files)]

# -----------------------------------------------------
# SPECIFICATIONS FOR THE USED COSMOLOGY
# -----------------------------------------------------
omega_b=0.045
omega_m=0.3
omega_l_in=0.7
pk_norm=0.81 
n=0.961

# -----------------------------------------------------
# UCAT / UFIG PARAMETERS
# -----------------------------------------------------
lum_fct_z_res = 0.001
lum_fct_m_res = 0.001 

lum_fct_alpha_blue = -1.3            
lum_fct_alpha_red = -0.5

template_coeff_z1_blue = 1                                     # Redshift z1>0 for blue galaxies
template_coeff_z1_red = 1                                      # Redshift z1>0 for red galaxies
magnitude_calculation = 'table'                                # The way magnitudes are calculated
lum_fct_filter_band = 'GenericBessel_B'                        # Filter band in which the luminosity function is valid

maps_remote_dir = 'ufig_res/maps/'
filters_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/filters_collection.h5'
templates_int_tables_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/sed_integrals__template_spectra_BlantonRoweis07.h5'

# PARAMETERS FROM COMMON.PY
lum_fct_m_star_blue_slope = -0.9408582               # Parameter a in M*(z) = a*z + b for blue galaxies, M*: Schechter parameter
lum_fct_m_star_blue_intcpt = -20.40492365            # Parameter b in M*(z) = a*z + b for blue galaxies, M*: Schechter parameter
lum_fct_m_star_red_slope = -0.70798041               # Parameter a in M*(z) = a*z + b for red galaxies, M*: Schechter parameter
lum_fct_m_star_red_intcpt = -20.37196157             # Parameter b in M*(z) = a*z + b for red galaxies, M*: Schechter parameter
lum_fct_phi_star_blue_amp = 0.00370253               # Parameter a in phi*(z) = a * exp(bz) for blue galaxies, phi*: Schechter parameter
lum_fct_phi_star_blue_exp = -0.10268436              # Parameter b in phi*(z) = a * exp(bz) for blue galaxies, phi*: Schechter parameter
lum_fct_phi_star_red_amp = 0.0035097                 # Parameter a in phi*(z) = a * exp(bz) for red galaxies, phi*: Schechter parameter
lum_fct_phi_star_red_exp = -0.70596888               # Parameter b in phi*(z) = a * exp(bz) for red galaxies, phi*: Schechter parameter

template_coeff_alpha0_blue_0 = 1.9946549                       # Dirichlet parameter for blue galaxies at z=0
template_coeff_alpha0_blue_1 = 1.99469164                      # Dirichlet parameter for blue galaxies at z=0
template_coeff_alpha0_blue_2 = 1.99461187                      # Dirichlet parameter for blue galaxies at z=0
template_coeff_alpha0_blue_3 = 1.9946589                       # Dirichlet parameter for blue galaxies at z=0
template_coeff_alpha0_blue_4 = 1.99463069                      # Dirichlet parameter for blue galaxies at z=0
template_coeff_alpha1_blue_0 = template_coeff_alpha0_blue_0    # Dirichlet parameter for blue galaxies at z=z1
template_coeff_alpha1_blue_1 = template_coeff_alpha0_blue_1    # Dirichlet parameter for blue galaxies at z=z1
template_coeff_alpha1_blue_2 = template_coeff_alpha0_blue_2    # Dirichlet parameter for blue galaxies at z=z1
template_coeff_alpha1_blue_3 = template_coeff_alpha0_blue_3    # Dirichlet parameter for blue galaxies at z=z1
template_coeff_alpha1_blue_4 = template_coeff_alpha0_blue_4    # Dirichlet parameter for blue galaxies at z=z1

template_coeff_alpha0_red_0 = 1.62158197                       # Dirichlet parameter for red galaxies at z=0
template_coeff_alpha0_red_1 = 1.62137391                       # Dirichlet parameter for red galaxies at z=0
template_coeff_alpha0_red_2 = 1.62175061                       # Dirichlet parameter for red galaxies at z=0
template_coeff_alpha0_red_3 = 1.62159144                       # Dirichlet parameter for red galaxies at z=0
template_coeff_alpha0_red_4 = 1.62165971                       # Dirichlet parameter for red galaxies at z=0
template_coeff_alpha1_red_0 = template_coeff_alpha0_red_0      # Dirichlet parameter for red galaxies at z=z1
template_coeff_alpha1_red_1 = template_coeff_alpha0_red_1      # Dirichlet parameter for red galaxies at z=z1
template_coeff_alpha1_red_2 = template_coeff_alpha0_red_2      # Dirichlet parameter for red galaxies at z=z1
template_coeff_alpha1_red_3 = template_coeff_alpha0_red_3      # Dirichlet parameter for red galaxies at z=z1
template_coeff_alpha1_red_4 = template_coeff_alpha0_red_4      # Dirichlet parameter for red galaxies at z=z1

# Weights for blue and red galaxies applied after drawing the coefficients
template_coeff_weight_blue = np.array([3.47116583e+09, 3.31262983e+06, 2.13298069e+09, 1.63722853e+10, 1.01368664e+09])
template_coeff_weight_red = np.array([3.84729278e+09, 1.56768931e+06, 3.91242928e+08, 4.66363319e+10, 3.03275998e+07])
scale_5500A = np.array([7.9213710e-12, 5.0803934e-09, 8.0514914e-11, 1.0282473e-11, 1.3742084e-10])
template_coeff_weight_blue *= scale_5500A
template_coeff_weight_red *= scale_5500A

m_star_par_blue = (lum_fct_m_star_blue_slope, 
                   lum_fct_m_star_blue_intcpt)
m_star_par_red = (lum_fct_m_star_red_slope,
                  lum_fct_m_star_red_intcpt)
phi_star_par_blue = (lum_fct_phi_star_blue_amp,
                     lum_fct_phi_star_blue_exp)
phi_star_par_red = (lum_fct_phi_star_red_amp,
                    lum_fct_phi_star_red_exp)

# Template coefficient parameters
template_coeff_alpha0_blue = np.array([template_coeff_alpha0_blue_0,
                                       template_coeff_alpha0_blue_1,
                                       template_coeff_alpha0_blue_2,
                                       template_coeff_alpha0_blue_3,
                                       template_coeff_alpha0_blue_4])
template_coeff_alpha1_blue = np.array([template_coeff_alpha1_blue_0,
                                       template_coeff_alpha1_blue_1,
                                       template_coeff_alpha1_blue_2,
                                       template_coeff_alpha1_blue_3,
                                       template_coeff_alpha1_blue_4])
template_coeff_alpha0_red = np.array([template_coeff_alpha0_red_0,
                                      template_coeff_alpha0_red_1,
                                      template_coeff_alpha0_red_2,
                                      template_coeff_alpha0_red_3,
                                      template_coeff_alpha0_red_4])
template_coeff_alpha1_red = np.array([template_coeff_alpha1_red_0,
                                      template_coeff_alpha1_red_1,
                                      template_coeff_alpha1_red_2,
                                      template_coeff_alpha1_red_3,
                                      template_coeff_alpha1_red_4])

# ----------------------------------------------------
# SET PYCOSMO COSMOLOGY
# -----------------------------------------------------
cosmo = PyCosmo.Cosmo()
cosmo.set(omega_b=omega_b, omega_m=omega_m, omega_l_in=omega_l_in, pk_norm=pk_norm, n=n_param)
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

MAGNITUDES_CALCULATOR = {'table': MagCalculatorTable}
loop_filter_names = [lum_fct_filter_band, 'DECam_r', 'DECam_g', 'DECam_z', 'DECam_i', 'DECam_Y']
filter_wavelengths = io_util.load_from_hdf5(filters_file_name, [loop_filter_names[i]+'/lam' for i in range(len(loop_filter_names))], root_path=maps_remote_dir)
filter_amplitudes = io_util.load_from_hdf5(filters_file_name, [loop_filter_names[i]+'/amp' for i in range(len(loop_filter_names))], root_path=maps_remote_dir)
filter_lam = {loop_filter_names[i]: filter_wavelengths[i] for i in range(len(loop_filter_names))}
filter_amp = {loop_filter_names[i]: filter_amplitudes[i] for i in range(len(loop_filter_names))}

loop_filters = Filter(lam=filter_lam, amp=filter_amp)
mag_calc = MAGNITUDES_CALCULATOR[magnitude_calculation](loop_filters)
n_templates = mag_calc.n_templates_dict[lum_fct_filter_band]

# ----------------------------------------------------–
# LOAD GALAXIES
# -----------------------------------------------------
abs_mag = np.load(dirname + infile_ucat_absmag)
z_ucat = np.load(dirname + infile_ucat_z)
blue_red = np.load(dirname + infile_ucat_redblue)

abs_mag_red = abs_mag[blue_red < 0.5]
z_red = z_ucat[blue_red < 0.5]
abs_mag_blue = abs_mag[blue_red > 0.5]
z_blue = z_ucat[blue_red > 0.5]

# ----------------------------------------------------–
# LOAD HALO-SUBHALO HISTOGRAMS
# -----------------------------------------------------
with load(infile_halos_dir + infile_hist_halos) as data:
	hist_z_mass_halos=data['hist_z_mass_halos']
	bin_edges_z=data['bin_edges_z']
	bin_edges_mass=data['bin_edges_mass']

with load(infile_halos_dir + infile_hist_subhalos) as data:
	hist_z_mass_subs=data['hist_z_mass_subs']

num_z_bins = len(bin_edges_z) - 1
num_mass_bins = len(bin_edges_mass) - 1

# ----------------------------------------------------–
# CREATE 2D HISTOGRAMS FOR RED vs. BLUE
# -----------------------------------------------------
# find index to match M_limit to a value in bin_edges_mass
idx_lim = (np.abs(bin_edges_mass - M_limit)).argmin()
M_limit_effective = bin_edges_mass[idx_lim]

mask_hist_red = np.concatenate((np.ones(idx_lim+1), np.zeros(num_mass_bins-idx_lim)))
mask_hist_blue = np.concatenate((np.zeros(idx_lim+1), np.ones(num_mass_bins-idx_lim)))

hist_z_mass_blue = hist_z_mass_halos * mask_hist_blue
hist_z_mass_red = hist_z_mass_halos * mask_hist_red + hist_z_mass_subs

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
np.savez(dirname + output_interp_red, lim_abs_mag_red=lim_abs_mag_red, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)
np.savez(dirname + output_interp_blue, lim_abs_mag_blue=lim_abs_mag_blue, bin_edges_z=bin_edges_z, bin_edges_mass=bin_edges_mass)

# ----------------------------------------------------–
# PREPARE STACKED MASS EDGES
# -----------------------------------------------------
mass_edges_stacked = np.reshape(np.tile(bin_edges_mass, len(bin_edges_z)), (len(bin_edges_mass)*len(bin_edges_z),))
z_edges_stacked = np.reshape(np.ndarray.flatten(np.transpose(np.tile(bin_edges_z, (len(bin_edges_mass), 1)))), (len(bin_edges_mass)*len(bin_edges_z),))

# ----------------------------------------------------–
# PREPARE EMPTY ARRAYS FOR GALAXIES
# -----------------------------------------------------
z = np.array([])
abs_mag = np.array([])
blue_red = np.array([])
app_mag_g = np.array([])
app_mag_r = np.array([])
app_mag_i = np.array([])
app_mag_z = np.array([])
app_mag_Y = np.array([])
halo_mass = np.array([])
x_coord = np.array([])
y_coord = np.array([])
z_coord = np.array([])
host_sub_index = np.array([])

n_uncut = 0
n_blue_uncut = 0
# ----------------------------------------------------–
# LOOP OVER HALO-SUBHALO FILES
# -----------------------------------------------------
for i in range(num_files):
	# ----------------------------------------------------–
	# LOAD HALO-SUBHALO FILE
	# -----------------------------------------------------
	mass, z_pin, x_coord_pin, y_coord_pin, z_coord_pin, host_sub = np.loadtxt(infile_halos_dir + infile_halos + file_ending[i] + '.txt', usecols=(1,2,3,4,5,-1), unpack=True)
	# ----------------------------------------------------–
	# DEVIDE HALOS AND SUBHALOS INTO RED AND BLUE
	# -----------------------------------------------------
	mask_red = (host_sub < 0.5) | (mass > M_limit_effective)
	mask_blue = ~ mask_red
	n_blue = len(z_pin[mask_blue])

	n_blue_uncut += n_blue
	n_uncut += len(z_pin)
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
	blue_red_temp = np.ones_like(z, dtype=np.int16)
	blue_red_temp[n_blue:] = 0
	halo_mass_temp = np.append(mass[mask_blue], mass[mask_red])
	x_coord_temp = np.append(x_coord_pin[mask_blue], x_coord_pin[mask_red])
	y_coord_temp = np.append(y_coord_pin[mask_blue], y_coord_pin[mask_red])
	z_coord_temp = np.append(z_coord_pin[mask_blue], z_coord_pin[mask_red])
	host_sub_index_temp = np.append(host_sub[mask_blue], host_sub[mask_red])
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
	# I'm setting them to zero, since I'm ignoring extionction values, as I don't have pixel coordinates
	excess_b_v = np.zeros(len(z_temp))
	# ----------------------------------------------------–
	# CALCULATE APPARENT MAGNITUDES
	# -----------------------------------------------------
	app_mag_g_temp = mag_calc(z_temp, excess_b_v, template_coeffs, ['DECam_g'])['DECam_g']
	app_mag_r_temp = mag_calc(z_temp, excess_b_v, template_coeffs, ['DECam_r'])['DECam_r']
	app_mag_i_temp = mag_calc(z_temp, excess_b_v, template_coeffs, ['DECam_i'])['DECam_i']
	app_mag_z_temp = mag_calc(z_temp, excess_b_v, template_coeffs, ['DECam_z'])['DECam_z']
	app_mag_Y_temp = mag_calc(z_temp, excess_b_v, template_coeffs, ['DECam_Y'])['DECam_Y']
	# ----------------------------------------------------–
	# APPLY CUTS IN APPARENT MAGNITUDES
	# -----------------------------------------------------
	mask_mag_range = (app_mag_i_temp >= gals_mag_min) & (app_mag_i_temp <= gals_mag_max)
	z_temp = z_temp[mask_mag_range]
	abs_mag_temp = abs_mag_temp[mask_mag_range]
	blue_red_temp = blue_red_temp[mask_mag_range]
	app_mag_g_temp = app_mag_g_temp[mask_mag_range]
	app_mag_r_temp = app_mag_r_temp[mask_mag_range]
	app_mag_i_temp = app_mag_i_temp[mask_mag_range]
	app_mag_z_temp = app_mag_z_temp[mask_mag_range]
	app_mag_Y_temp = app_mag_Y_temp[mask_mag_range]
	halo_mass_temp = halo_mass_temp[mask_mag_range]
	x_coord_temp = x_coord_temp[mask_mag_range]
	y_coord_temp = y_coord_temp[mask_mag_range]
	z_coord_temp = z_coord_temp[mask_mag_range]
	host_sub_index_temp = host_sub_index_temp[mask_mag_range]
	# ----------------------------------------------------–
	# APPEND REMAINING GALAXIES
	# -----------------------------------------------------
	z = np.append(z, z_temp)
	abs_mag = np.append(abs_mag, abs_mag_temp)
	blue_red = np.append(blue_red, blue_red_temp)
	app_mag_g = np.append(app_mag_g, app_mag_g_temp)
	app_mag_r = np.append(app_mag_r, app_mag_r_temp)
	app_mag_i = np.append(app_mag_i, app_mag_i_temp)
	app_mag_z = np.append(app_mag_z, app_mag_z_temp)
	app_mag_Y = np.append(app_mag_Y, app_mag_Y_temp)
	halo_mass = np.append(halo_mass, halo_mass_temp)
	x_coord = np.append(x_coord, x_coord_temp)
	y_coord = np.append(y_coord, y_coord_temp)
	z_coord = np.append(z_coord, z_coord_temp)
	host_sub_index = np.append(host_sub_index, host_sub_index_temp)
# end of loop
# ----------------------------------------------------–
# SAVE OUTPUT
# -----------------------------------------------------
# each property in a separate npy file
np.save(dirname+outfile_galaxies_base+'z', z)
np.save(dirname+outfile_galaxies_base+'abs_mag', abs_mag)
np.save(dirname+outfile_galaxies_base+'blue_red', blue_red)
np.save(dirname+outfile_galaxies_base+'app_mag_g', app_mag_g)
np.save(dirname+outfile_galaxies_base+'app_mag_r', app_mag_r)
np.save(dirname+outfile_galaxies_base+'app_mag_i', app_mag_i)
np.save(dirname+outfile_galaxies_base+'app_mag_z', app_mag_z)
np.save(dirname+outfile_galaxies_base+'app_mag_Y', app_mag_Y)
np.save(dirname+outfile_galaxies_base+'halo_mass', halo_mass)
np.save(dirname+outfile_galaxies_base+'x_coord', x_coord)
np.save(dirname+outfile_galaxies_base+'y_coord', y_coord)
np.save(dirname+outfile_galaxies_base+'z_coord', z_coord)
np.save(dirname+outfile_galaxies_base+'host_sub_index', host_sub_index)

# ----------------------------------------------------–
# SOME PRINT STATEMENTS
# -----------------------------------------------------

print('Code run.')
print('num_files = ' + str(num_files))
print('infile_halos_dir = ' + infile_halos_dir)
print('infile_hist_halos = ' + infile_hist_halos)
print('infile_hist_subhalos = ' + infile_hist_subhalos)
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
print('omega_l_in = ' + str(omega_l_in))
print('pk_norm = ' + str(pk_norm))
print('n = ' + str(n))
print('filters = ' + str([loop_filter_names[i] for i in range(len(loop_filter_names))]))
print('total number of galaxies without magnitude cuts = ' + str(n_uncut))
print('number of blue galaxies without magnitude cuts = ' + str(n_blue_uncut))
print('number of red galaxies without magnitude cuts = ' + str(n_uncut - n_blue_uncut))
print('total number of galaxies after magnitude cuts = ' + str(len(z)))
print('number of blue galaxies after magnitude cuts = ' + str(np.sum(host_sub_index)))
print('number of red galaxies after magnitude cuts = ' + str(len(z) - np.sum(host_sub_index)))
