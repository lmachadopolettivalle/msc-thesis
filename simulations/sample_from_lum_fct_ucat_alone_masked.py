# Description of the program:
# sample galaxies
# calculate apparent magnitudes, apply cuts, save output
# Note: with a spatial mask

# Author: Pascale Berner
# first written: 15.12.2022
# last adapted: 15.12.2022
# partially copied from: sample_from_lum_fct_only_masked.py and sample_from_lum_fct_ucat_alone_sdss_joerg.py

# ----------------------------------------------------
# IMPORTS
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
import pandas as pd
import healpy as hp

# ----------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------

z_max = 0.75
m_max = -12

DEFAULT_PIXAREA = 0.1*4*np.pi # set pixarea to the desired value, if no mask is provided, otherwise leave it
#DEFAULT_PIXAREA = 0.00001*4*np.pi
infile_footprint = '/cluster/work/refregier/bernerp/DES/DES_Y1/y1a1_gold_1.0.2_wide_footprint_4096.fit'
NEST=False
MASK_CUTOFF = 0
dec_min = -90
dec_max = -30

dirname = './'
outfile_galaxies_base = 'sampled_app_mag_desy1_20221215_1_'

# apparent magnitude limits, in the main band (typically i or r)
gals_mag_max = 20.
gals_mag_min = 10.

# -----------------------------------------------------
# SPECIFICATIONS FOR THE USED COSMOLOGY
# -----------------------------------------------------
omega_b=0.045
omega_m=0.276
omega_l_in=0.724
pk_norm=0.811 
n_param=0.961

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

print('params defined')

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

# ----------------------------------------------------
# LOAD MASK, CALCULATE PIXAREA
# -----------------------------------------------------
pixarea = DEFAULT_PIXAREA

m_footprint = hp.fitsfunc.read_map(infile_footprint)
if m_footprint is not None:
	NSIDE = hp.get_nside(m_footprint)
	NPIX = hp.nside2npix(NSIDE)

	higher_colatitude = np.pi/2 - np.radians(dec_min)
	lower_colatitude = np.pi/2 - np.radians(dec_max)

	ring_dec_range_pixel_ids = hp.query_strip(NSIDE, theta1=lower_colatitude, theta2=higher_colatitude, inclusive=False, nest=False)
	dec_range_pixel_ids = ring_dec_range_pixel_ids

	if NEST is True:
		ring_dec_range_map = np.zeros(NPIX)
		ring_dec_range_map[ring_dec_range_pixel_ids] = 1
		nest_dec_range_map = hp.reorder(ring_dec_range_map, r2n=True)
		nest_dec_range_pixel_ids = np.where(nest_dec_range_map > 0)[0]
		dec_range_pixel_ids = nest_dec_range_pixel_ids

	mask_pixel_ids = np.where(m_footprint > MASK_CUTOFF)[0]
	individual_pixarea = hp.pixelfunc.nside2pixarea(NSIDE, degrees=False)
	filtered_pixel_ids = np.intersect1d(mask_pixel_ids, dec_range_pixel_ids)

	pixarea = individual_pixarea * len(filtered_pixel_ids)

print('mask loaded')

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
index_blue = np.ones(len(z_blue)) # 1 = blue
index_red = np.zeros(len(z_red)) # 0 = red

n_tot = n_gal_red + n_gal_blue
z_temp = np.concatenate((z_blue, z_red))
abs_mag_temp = np.concatenate((abs_mag_blue, abs_mag_red))
blue_red_temp = np.concatenate((index_blue, index_red))

n_blue = len(z_blue)

print('galaxies sampled')

# ----------------------------------------------------–
# CALCULATE TEMPLATE COEFFICIENTS
# -----------------------------------------------------
# Draw template coefficients
template_coeffs = np.empty((len(z_temp), n_templates))
template_coeffs[:n_gal_blue] = galaxy_sed.sample_template_coeff_dirichlet(z_blue, template_coeff_alpha0_blue, template_coeff_alpha1_blue, template_coeff_z1_blue, template_coeff_weight_blue)
template_coeffs[n_gal_blue:] = galaxy_sed.sample_template_coeff_dirichlet(z_red, template_coeff_alpha0_red, template_coeff_alpha1_red, template_coeff_z1_red, template_coeff_weight_red)
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

# ----------------------------------------------------–
# SAVE OUTPUT
# -----------------------------------------------------
# each property in a separate npy file
np.save(dirname+outfile_galaxies_base+'z', z_temp)
np.save(dirname+outfile_galaxies_base+'abs_mag', abs_mag_temp)
np.save(dirname+outfile_galaxies_base+'blue_red', blue_red_temp)
np.save(dirname+outfile_galaxies_base+'app_mag_g', app_mag_g_temp)
np.save(dirname+outfile_galaxies_base+'app_mag_r', app_mag_r_temp)
np.save(dirname+outfile_galaxies_base+'app_mag_i', app_mag_i_temp)
np.save(dirname+outfile_galaxies_base+'app_mag_z', app_mag_z_temp)
np.save(dirname+outfile_galaxies_base+'app_mag_Y', app_mag_Y_temp)

# ----------------------------------------------------–
# SOME PRINT STATEMENTS
# -----------------------------------------------------

print('Code run.')

print('gals_mag_max = ' + str(gals_mag_max))
print('gals_mag_min = ' + str(gals_mag_min))
print('omega_b = ' + str(omega_b))
print('omega_m = ' + str(omega_m))
print('omega_l_in = ' + str(omega_l_in))
print('pk_norm = ' + str(pk_norm))
print('n = ' + str(n_param))
print('filters = ' + str([loop_filter_names[i] for i in range(len(loop_filter_names))]))
print('total number of galaxies without magnitude cuts = ' + str(n_gal_blue + n_gal_red))
print('number of blue galaxies without magnitude cuts = ' + str(n_gal_blue))
print('number of red galaxies without magnitude cuts = ' + str(n_gal_red))
print('total number of galaxies after magnitude cuts = ' + str(len(z_temp)))
print('number of blue galaxies after magnitude cuts = ' + str(np.sum(blue_red_temp)))
print('number of red galaxies after magnitude cuts = ' + str(len(z_temp) - np.sum(blue_red_temp)))




