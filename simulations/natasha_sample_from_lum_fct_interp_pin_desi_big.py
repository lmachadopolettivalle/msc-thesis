# LUIS: tested with modifications on 17.11


# Sample z and M_B from luminosity functions
# get z and M_halo from PINOCCHIO (including subhalos) in lightcone
# 2D interpolation to get M_B(z, M_halo)
# get apparent magnitudes

# required inputs: halo-subhalo catalog in lightcone, matching galaxy catalog

# author: Pascale Berner
# author2: Natasha Theiler
# last changed: 16.05.2021



# ---------------------------------------
# IMPORTS

import time
start_time = time.clock()

from ucat import galaxy_sampling_util
from ucat import io_util
from ucat.galaxy_population_models import galaxy_luminosity_function
from ucat.galaxy_population_models import galaxy_sed
#from ucat.plugins import sample_galaxies
from scipy.interpolate import griddata
import h5py
import scipy.integrate
import PyCosmo
import numpy as np

print('imports done')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# DEFINE CLASSES

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
        self.z_grid, self.excess_b_v_grid = io_util.load_from_hdf5(templates_int_tables_file_name,
                                                                   ('z', 'E(B-V)'),
                                                                   root_path=maps_remote_dir)
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

#            mags = -2.5 * np.log10(mags / (self.c * self.filter_norm_dict[filter_name])) - 48.6
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
                                        range(self.n_templates_dict[filter_name])],
                                       root_path=self.maps_remote_dir)
            self.filter_norm_dict[filter_name] = scipy.integrate.simps(filters.amp[filter_name] / filters.lam[filter_name],
                                                                       x=filters.lam[filter_name])

    def get_n_templates(self, filter_names):
        with h5py.File(io_util.get_abs_path(self.templates_int_tables_file_name, root_path=self.maps_remote_dir),
                       mode='r') as f:
            self.n_templates_dict = {filter_name: len(list(filter(lambda k: k.startswith('template_'),
                                                                  f['integrals/{}'.format(filter_name)].keys())))
                                     for filter_name in filter_names}

print('classes defined')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# INPUTS
# ADAPT THIS SECTION!

# NOTE: AREAS OF LIGHCONE AND UCAT FILE HAVE TO MATCH. IF DEPTH DOES NOT MATCH: APPLY A CUT AFTERWARDS.

#pinocchio_halo_subhalo_catalog = '/cluster/scratch/bernerp/pinocchio/pin_500mpc_2048_9/pin_500mpc_2048_9_halo_subhalo_sorted_plc_new_4_1.txt'
#pinocchio_halo_subhalo_catalog = '/cluster/work/refregier/bernerp/DESI_BGS_filters/pin_500mpc_2048_9_halo_subhalo_sorted_plc_new_4_1_zlt05.txt'
#pinocchio_halo_subhalo_catalog = '/cluster/work/refregier/bernerp/DESI_BGS_filters/pin_500mpc_2048_8_halo_subhalo_sorted_plc_new_6_1.txt'
#pinocchio_halo_subhalo_catalog = '/cluster/work/refregier/bernerp/DESI_BGS_filters/pin_500mpc_2048_9_halo_subhalo_sorted_plc_new_6_1.txt'
#pinocchio_halo_subhalo_catalog = '/cluster/work/refregier/bernerp/DESI_BGS_filters/pin_500mpc_2048_8_halo_subhalo_sorted_plc_new_6_2.txt'
pinocchio_halo_subhalo_catalog = '/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/pinocchio_masked_halos_subhalos_plc.0.txt'

#z_max = 1.5
z_max = 0.5
z_min = 0.0
m_max = -12.3

#infile_ucat = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_7.5_20220302.txt'
infile_ucat_z = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_60_z_tot_20220411_123.npy'

infile_ucat_abs_mag = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_60_abs_mag_tot_20220411_123.npy'

infile_ucat_index = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_60_index_tot_20220411_123.npy'

#infile_ucat_z = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_15_z_tot_20220430_123.npy'

#infile_ucat_abs_mag = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_15_abs_mag_tot_20220430_123.npy'

#infile_ucat_index = '/cluster/scratch/ntheiler/ucat/ucat_sorted_only_theta_15_index_tot_20220430_123.npy'

#outfile_galaxies = '/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220302_pin500mpc20489_zlt05.txt'
#outfile_galaxies = '/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220323_pin500mpc20489_zlt05_13_m.txt'

#M_limit = 8.e12 # mass limit for assigning blue or red galaxies to halos, [Msun/h]
M_limit = 8.0*10**13
num_z_bins = 150  # number of bins for interpolation
num_mass_bins = 30

omega_b=0.045
omega_m=0.3
omega_l_in=0.7
pk_norm=0.81 
n=0.961


# ---------------------------------------
# PARAMTERS
# ADAPT THIS SECTION!

lum_fct_z_res = 0.001
lum_fct_m_res = 0.001
lum_fct_alpha_blue = -1.3
lum_fct_alpha_red = -0.5

lum_fct_m_star_blue_slope = -0.9408582               # Parameter a in M*(z) = a*z + b for blue galaxies, M*: Schechter parameter
lum_fct_m_star_blue_intcpt = -20.40492365            # Parameter b in M*(z) = a*z + b for blue galaxies, M*: Schechter parameter
lum_fct_m_star_red_slope = -0.70798041               # Parameter a in M*(z) = a*z + b for red galaxies, M*: Schechter parameter
lum_fct_m_star_red_intcpt = -20.37196157             # Parameter b in M*(z) = a*z + b for red galaxies, M*: Schechter parameter
lum_fct_phi_star_blue_amp = 0.00370253               # Parameter a in phi*(z) = a * exp(bz) for blue galaxies, phi*: Schechter parameter
lum_fct_phi_star_blue_exp = -0.10268436              # Parameter b in phi*(z) = a * exp(bz) for blue galaxies, phi*: Schechter parameter
lum_fct_phi_star_red_amp = 0.0035097                 # Parameter a in phi*(z) = a * exp(bz) for red galaxies, phi*: Schechter parameter
lum_fct_phi_star_red_exp = -0.70596888               # Parameter b in phi*(z) = a * exp(bz) for red galaxies, phi*: Schechter parameter

template_coeff_z1_blue = 1                                     # Redshift z1>0 for blue galaxies
template_coeff_z1_red = 1                                      # Redshift z1>0 for red galaxies
magnitude_calculation = 'table'                                # The way magnitudes are calculated
lum_fct_filter_band = 'GenericBessel_B'                                      # Filter band in which the luminosity function is valid

#filters_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/template_spectra_BlantonRoweis07.h5'  # had to adapt
filters_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/filters_collection.h5'

maps_remote_dir = 'ufig_res/maps/'                             # have to adapt                          # Remote directory containing maps
templates_int_tables_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/sed_integrals__template_spectra_BlantonRoweis07.h5'  # had to adapt
#par_filters = ['g', 'r', 'i', 'z', 'Y']                        # Filter bands (multi-band only)

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

template_coeff_weight_blue = np.array([3.47116583e+09, 3.31262983e+06, 2.13298069e+09, 1.63722853e+10, 1.01368664e+09])
template_coeff_weight_red = np.array([3.84729278e+09, 1.56768931e+06, 3.91242928e+08, 4.66363319e+10, 3.03275998e+07])
scale_5500A = np.array([7.9213710e-12, 5.0803934e-09, 8.0514914e-11, 1.0282473e-11, 1.3742084e-10])
template_coeff_weight_blue *= scale_5500A
template_coeff_weight_red *= scale_5500A

#gals_mag_max = 28
gals_mag_max = 19.5
gals_mag_min = 15

# Luminosity function parameters
m_star_par_blue = (lum_fct_m_star_blue_slope, lum_fct_m_star_blue_intcpt)
m_star_par_red = (lum_fct_m_star_red_slope, lum_fct_m_star_red_intcpt)
phi_star_par_blue = (lum_fct_phi_star_blue_amp, lum_fct_phi_star_blue_exp)
phi_star_par_red = (lum_fct_phi_star_red_amp, lum_fct_phi_star_red_exp)

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

print('parameters set')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# COSMOLOGY

cosmo = PyCosmo.Cosmo()
cosmo.set(omega_b=omega_b, omega_m=omega_m, omega_l_in=omega_l_in, pk_norm=pk_norm, n=n)

print('Cosmology set')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# GET GALAXIES FROM A PRE-SAMPLED CATALOG
# using: sample_from_lum_fct_only.py

#abs_mag, z_ucat, blue_red = np.loadtxt(infile_ucat, usecols=(0,1,2), unpack=True)
"""
NOTE changed this
abs_mag = np.load(infile_ucat_abs_mag, allow_pickle=True)
z_ucat = np.load(infile_ucat_z, allow_pickle=True)
blue_red = np.load(infile_ucat_index, allow_pickle=True)
"""
blue_red = np.random.choice(2, 1000)
abs_mag = -10 * np.ones(len(blue_red))
z_ucat = 0.4 * np.ones(len(blue_red))


abs_mag_red = abs_mag[blue_red < 0.5]
z_red = z_ucat[blue_red < 0.5]
abs_mag_blue = abs_mag[blue_red > 0.5]
z_blue = z_ucat[blue_red > 0.5]

del z_ucat
del abs_mag
del blue_red

print('Galaxy catalog loaded from file')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

# ---------------------------------------
# LOAD HALO FILE, 2D HISTOGRAM

# host_sub = 1 for halos, 0 for subhalos
mass, z_pin, x_coord_pin, y_coord_pin, z_coord_pin, host_sub = np.loadtxt(pinocchio_halo_subhalo_catalog, usecols=(1,2,3,4,5,-1), unpack=True)
# IF THE LIGHTCONE IS TOO DEEP: APPLY A CUT IN Z_PIN

print('halos loaded')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

# mask to assign blue or red galaxies to the halos/subhalos
mask_red = (host_sub < 0.5) | (mass > M_limit)
mask_blue = ~ mask_red


# define bins for 2D histogram
min_mass = min(mass)
max_mass = max(mass)


bin_edges_mass = np.logspace(np.log10(min_mass), np.log10(max_mass), num=(num_mass_bins+1))
print(bin_edges_mass)
bin_edges_z = np.linspace(z_min, z_max, (num_z_bins+1))

# 2D histogram for halos/subhalos (for red/blue separately)
z_pin_red = z_pin[mask_red]
z_pin_blue = z_pin[mask_blue]

del z_pin

hist_z_mass_red, bin_edges_z, bin_edges_mass = np.histogram2d(z_pin_red, mass[mask_red], bins=(bin_edges_z, bin_edges_mass))
hist_z_mass_blue, bin_edges_z, bin_edges_mass = np.histogram2d(z_pin_blue, mass[mask_blue], bins=(bin_edges_z, bin_edges_mass))

print('2D hist done')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# LOOPS TO GET VALUES FOR INTERPOLATION

# empty arrays for values; limit meaning at the edges of 2D histogram
lim_abs_mag_red = np.zeros((num_z_bins+1, num_mass_bins+1))
lim_abs_mag_blue = np.zeros((num_z_bins+1, num_mass_bins+1))
lim_abs_mag_red[:] = np.nan
lim_abs_mag_blue[:] = np.nan

num_z = num_z_bins + 1
num_mass = num_mass_bins + 1

array_leftover_blue = np.zeros(num_z_bins)
array_leftover_red = np.zeros(num_z_bins)

for i in range(num_z - 1):  # loop over redshift, starting at low z
	abs_mag_red_i = abs_mag_red[(z_red > bin_edges_z[i]) & (z_red <= bin_edges_z[i+1])]
	abs_mag_blue_i = abs_mag_blue[(z_blue > bin_edges_z[i]) & (z_blue <= bin_edges_z[i+1])]

	ind_red = 0
	ind_blue = 0

	for j in reversed(range(num_mass - 1)):  # loop over mass bins, starting at hight mass
		num_blue = hist_z_mass_blue[i,j]
		if int(ind_blue + num_blue + 1) >= len(abs_mag_blue_i):  # ensure there are enough galaxies
                        ind_blue = int(ind_blue + num_blue)
                        break
		if num_blue > 0:
			#mean_abs_mag_blue[i,j] = np.mean(abs_mag_blue_i[ind_blue : int(ind_blue + num_blue)])
			lim_abs_mag_blue[i+1,j] = abs_mag_blue_i[int(ind_blue + num_blue)]
			if np.isnan(lim_abs_mag_blue[i,j]):  # check if the neighboring bin in redshift is empty
				lim_abs_mag_blue[i,j] = abs_mag_blue_i[int(ind_blue + num_blue)]
			if np.isnan(lim_abs_mag_blue[i+1,j+1]):  # check if the neighboring bin in mass is empty
				lim_abs_mag_blue[i+1,j+1] = abs_mag_blue_i[ind_blue]
			ind_blue =  int(ind_blue + num_blue)
	array_leftover_blue[i] = (len(abs_mag_blue_i) - ind_blue)

	for j in reversed(range(num_mass - 1)):  # loop over mass bins, starting at hight mass
		num_red = hist_z_mass_red[i,j]
		if int(ind_red + num_red + 1) >= len(abs_mag_red_i):  # ensure there are enough galaxies
			ind_red = int(ind_red + num_red)
			break
		if num_red > 0:
			#mean_abs_mag_red[i,j] = np.mean(abs_mag_red_i[ind_red : int(ind_red + num_red)])
			lim_abs_mag_red[i+1,j] = abs_mag_red_i[int(ind_red + num_red)]
			if np.isnan(lim_abs_mag_red[i,j]):  # check if the neighboring bin in redshift is empty
				lim_abs_mag_red[i,j] = abs_mag_red_i[int(ind_red + num_red)]
			if np.isnan(lim_abs_mag_red[i+1,j+1]):  # check if the neighboring bin in mass is empty
				lim_abs_mag_red[i+1,j+1] = abs_mag_red_i[ind_red]
			ind_red = int(ind_red + num_red)
	array_leftover_red[i] = (len(abs_mag_red_i) - ind_red)


print('array_leftover_blue', array_leftover_blue)
"""
NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_hist_array_leftover_blue', array_leftover_blue)
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_hist_array_leftover_red', array_leftover_red)
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_hist_array_leftover_bins', bin_edges_z)
"""

print('loop done')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

del hist_z_mass_red
del hist_z_mass_blue
del num_z
del num_mass
del abs_mag_blue
del z_blue
del z_red
del abs_mag_red_i
del abs_mag_blue_i
del num_blue
del num_red


#---------------------------------------
# 2D INTERPOLATION

mass_edges_stacked = np.reshape(np.tile(bin_edges_mass, len(bin_edges_z)), (len(bin_edges_mass)*len(bin_edges_z),))
z_edges_stacked = np.reshape(np.ndarray.flatten(np.transpose(np.tile(bin_edges_z, (len(bin_edges_mass), 1)))), (len(bin_edges_mass)*len(bin_edges_z),))


print('mass_edges_stacked nan: ', np.count_nonzero(np.isnan(mass_edges_stacked)))
print('z_edges_stacked nan: ', np.count_nonzero(np.isnan(z_edges_stacked)))
print('lim_abs_mag_red nan: ', np.count_nonzero(np.isnan(lim_abs_mag_red)))
print('z_pin_red nan: ', np.count_nonzero(np.isnan(z_pin_red)))
print('mass[mask_red]) nan: ', np.count_nonzero(np.isnan(mass[mask_red])))
print('lim_abs_mag_blue nan: ', np.count_nonzero(np.isnan(lim_abs_mag_blue)))
print('z_pin_blue nan: ', np.count_nonzero(np.isnan(z_pin_blue)))
print('mass[mask_blue]) nan: ', np.count_nonzero(np.isnan(mass[mask_blue])))

del bin_edges_z
del bin_edges_mass

abs_mag_red_pin = griddata((z_edges_stacked, np.log10(mass_edges_stacked)), np.ndarray.flatten(lim_abs_mag_red), (z_pin_red, np.log10(mass[mask_red])))
abs_mag_blue_pin = griddata((z_edges_stacked, np.log10(mass_edges_stacked)), np.ndarray.flatten(lim_abs_mag_blue), (z_pin_blue, np.log10(mass[mask_blue])))
del mass_edges_stacked
del z_edges_stacked
del lim_abs_mag_red
del lim_abs_mag_blue

z = np.append(z_pin_blue, z_pin_red)
abs_mag = np.append(abs_mag_blue_pin, abs_mag_red_pin)

z_hist_nans_blue, bin_edges_z_blue = np.histogram(z_pin_blue[np.isnan(abs_mag_blue_pin)], bins=50)
""" NOTE changed here
np.savez('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_z_hist_nans_blue', z_hist_nans_blue, bin_edges_z_blue)
"""

print('bin_edges_z_blue: ', bin_edges_z_blue)
print('z_pin_blue: ', z_pin_blue)
print('np.isnan(abs_mag_blue_pin): ', np.isnan(abs_mag_blue_pin))
""" NOTE changed here
np.savez('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_z_hist_nans_blue2', z_pin_blue, abs_mag_blue_pin)
"""


z_hist_nans_red, bin_edges_z_red = np.histogram(z_pin_red[np.isnan(abs_mag_red_pin)], bins=50)
""" NOTE changed here
np.savez('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220519_pin500mpc20489_4_1_zlt05_123_m_z_hist_nans_red', z_hist_nans_red, bin_edges_z_red)
"""


del abs_mag_red_pin
del abs_mag_blue_pin

# Create vector indicating which galaxy is blue and which is red: 1 = blue, 0 = red
blue_red = np.ones_like(z, dtype=np.int16)
blue_red[len(z_pin_blue):] = 0

print('interpolation done')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

print(len(blue_red))
print('Number of red galaxies with non-nan magnitudes: ', np.count_nonzero(~np.isnan(abs_mag) & (blue_red < 0.5)))
print('Number of blue galaxies with non-nan magnitudes: ', np.count_nonzero(~np.isnan(abs_mag) & (blue_red > 0.5)))
print('Number of red galaxies with nan magnitudes: ', np.count_nonzero(np.isnan(abs_mag) & (blue_red < 0.5)))
print('Number of blue galaxies with nan magnitudes: ', np.count_nonzero(np.isnan(abs_mag) & (blue_red > 0.5)))
print(len(abs_mag))



# ---------------------------------------
# DEFINE MAGNITUDE CALCULATOR

MAGNITUDES_CALCULATOR = {'table': MagCalculatorTable}
loop_filter_names = [lum_fct_filter_band, 'BASSMzLS_r', 'BASSMzLS_g', 'BASSMzLS_z', 'DECam_r', 'DECam_g', 'DECam_z']
filter_wavelengths = io_util.load_from_hdf5(filters_file_name, [loop_filter_names[i]+'/lam' for i in range(len(loop_filter_names))], root_path=maps_remote_dir)
filter_amplitudes = io_util.load_from_hdf5(filters_file_name, [loop_filter_names[i]+'/amp' for i in range(len(loop_filter_names))], root_path=maps_remote_dir)
filter_lam = {loop_filter_names[i]: filter_wavelengths[i] for i in range(len(loop_filter_names))}
filter_amp = {loop_filter_names[i]: filter_amplitudes[i] for i in range(len(loop_filter_names))}

loop_filters = Filter(lam=filter_lam, amp=filter_amp)
mag_calc = MAGNITUDES_CALCULATOR[magnitude_calculation](loop_filters)
n_templates = mag_calc.n_templates_dict[lum_fct_filter_band]


# ---------------------------------------
# TEMPLATE COEFFICIENTS

# Draw template coefficients
template_coeffs = np.empty((len(z), n_templates))
template_coeffs[:len(z_pin_blue)] = galaxy_sed.sample_template_coeff_dirichlet(z_pin_blue, template_coeff_alpha0_blue, template_coeff_alpha1_blue, template_coeff_z1_blue, template_coeff_weight_blue)
template_coeffs[len(z_pin_blue):] = galaxy_sed.sample_template_coeff_dirichlet(z_pin_red, template_coeff_alpha0_red, template_coeff_alpha1_red, template_coeff_z1_red, template_coeff_weight_red)

del z_pin_blue
del z_pin_red

# Calculate absolute magnitudes according to coefficients and adjust coefficients according to drawn magnitudes
template_coeffs *= np.expand_dims(10 ** (0.4 * (mag_calc(np.zeros_like(z),
                                                         np.zeros_like(z),
                                                         template_coeffs,
                                                         [lum_fct_filter_band])[lum_fct_filter_band]
                                                - abs_mag)), -1)

# Transform to apparent coefficients
lum_dist = galaxy_sampling_util.apply_pycosmo_distfun(cosmo.background.dist_lum_a, z)
template_coeffs *= np.expand_dims((10e-6 / lum_dist) ** 2 / (1 + z), -1)

del lum_dist

# ---------------------------------------
# EXTINCTION MAP

# exctinction_eval = sample_galaxies.ExtinctionMapEvaluator(par)
# excess_b_v = exctinction_eval(w, x, y)

excess_b_v = np.zeros(len(z))  # since I don't have pixel coordinates

print('prep for app mag done')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# ASSIGN HALO PROPERTIES (POSITIONS etc.)

halo_mass = np.append(mass[mask_blue], mass[mask_red])
x_coord = np.append(x_coord_pin[mask_blue], x_coord_pin[mask_red])
y_coord = np.append(y_coord_pin[mask_blue], y_coord_pin[mask_red])
z_coord = np.append(z_coord_pin[mask_blue], z_coord_pin[mask_red])
host_sub_index = np.append(host_sub[mask_blue], host_sub[mask_red])

del mask_blue
del mask_red
del x_coord_pin
del y_coord_pin
del z_coord_pin
del host_sub
del mass

# ---------------------------------------
# CALCULATE SPHERICAL COORDINATES

r_coord = np.sqrt(x_coord**2 + y_coord**2 + z_coord**2)
theta_coord = np.arctan2(np.sqrt(x_coord**2 + y_coord**2), z_coord)
#if x_coord == 0:
#   phi_coord = np.pi /2
#if x_coord != 0:
#   phi_coord = np.arctan2(y_coord, x_coord)
phi_coord = np.arctan2(y_coord, x_coord) + np.pi

if np.max(theta_coord > np.pi):
   print('theta wrong')
if np.max(phi_coord > np.pi):
   print('phi wrong')


# ---------------------------------------
# MASK FOR DECAM, BASS

#mask ratio
ratio_bass = 0.4655883621215696
bass_decam = phi_coord > 2*np.pi*ratio_bass
print('mask DECam:', bass_decam, 1*(bass_decam))

#mask lines
#angle1 = 3*np.pi/8
#angle2 = -3*np.pi/8
#bass_decam = (phi_coord < angle2) | (phi_coord > angle1)
#print('mask DECam:', bass_decam, 1*(bass_decam))

# ---------------------------------------
# CALCULATE APPARENT r-band MAGNITUDE

# note: use 'b' instead of 'r' if b-band is preferred
app_mag_r_b = mag_calc(z, excess_b_v, template_coeffs, ['BASSMzLS_r'])['BASSMzLS_r']
app_mag_g_b = mag_calc(z, excess_b_v, template_coeffs, ['BASSMzLS_g'])['BASSMzLS_g']
app_mag_z_b = mag_calc(z, excess_b_v, template_coeffs, ['BASSMzLS_z'])['BASSMzLS_z']
print('length appmag B: ', len(app_mag_r_b))

app_mag_r_d = mag_calc(z, excess_b_v, template_coeffs, ['DECam_r'])['DECam_r']
app_mag_g_d = mag_calc(z, excess_b_v, template_coeffs, ['DECam_g'])['DECam_g']
app_mag_z_d = mag_calc(z, excess_b_v, template_coeffs, ['DECam_z'])['DECam_z']
print('length appmag D: ', len(app_mag_r_d))

del excess_b_v
del template_coeffs

app_mag_r = (1-bass_decam)*app_mag_r_b + (bass_decam)*app_mag_r_d
del app_mag_r_b
del app_mag_r_d
app_mag_g = (1-bass_decam)*app_mag_g_b + (bass_decam)*app_mag_g_d
del app_mag_g_b
del app_mag_g_d
app_mag_z = (1-bass_decam)*app_mag_z_b + (bass_decam)*app_mag_z_d
del app_mag_z_b
del app_mag_z_d


#del app_mag_r_b
#del app_mag_g_b
#del app_mag_z_b
#del app_mag_r_d
#del app_mag_g_d
#del app_mag_z_d

print('length appmag r: ', len(app_mag_r))
print('length appmag z: ', len(app_mag_z))
print('length appmag g: ', len(app_mag_g))


print('app mag calculated')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# APPLY CUTS

print('Number of galaxies before cut: ' + str(len(abs_mag)))

mask_mag_range = (app_mag_r >= gals_mag_min) & (app_mag_r <= gals_mag_max)
n_gal = np.count_nonzero(mask_mag_range)
print('Number of galaxies after cut: ' + str(n_gal))
print('Cut at: app_mag < ' + str(gals_mag_max))
z = z[mask_mag_range]
abs_mag = abs_mag[mask_mag_range]
blue_red = blue_red[mask_mag_range]
app_mag_r = app_mag_r[mask_mag_range]
app_mag_g = app_mag_g[mask_mag_range]
app_mag_z = app_mag_z[mask_mag_range]

halo_mass = halo_mass[mask_mag_range]
x_coord = x_coord[mask_mag_range]
y_coord = y_coord[mask_mag_range]
z_coord = z_coord[mask_mag_range]
host_sub_index = host_sub_index[mask_mag_range]

r_coord = r_coord[mask_mag_range]
theta_coord = theta_coord[mask_mag_range]
phi_coord = phi_coord[mask_mag_range]
bass_decam = bass_decam[mask_mag_range]

del mask_mag_range

print('cuts applied')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

# GALAXIES WITH RANDOM POSITIONS AT LOW REDSHIFT TO FILL UP
# not yet done


# ---------------------------------------
# SAFE TO FILE

z = z.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_z', z)
"""
abs_mag = abs_mag.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_abs_mag', abs_mag)
"""
app_mag_g = app_mag_g.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_app_mag_g', app_mag_g)
"""
del app_mag_g
app_mag_r = app_mag_r.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_app_mag_r', app_mag_r)
"""
del app_mag_r
app_mag_z = app_mag_z.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_app_mag_z', app_mag_z)
"""
del app_mag_z
blue_red = blue_red.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_blue_red', blue_red)
"""
halo_mass = halo_mass.reshape((len(z),1))
x_coord = x_coord.reshape((len(z),1))
y_coord = y_coord.reshape((len(z),1))
z_coord = z_coord.reshape((len(z),1))
r_coord = r_coord.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_r_coord', r_coord)
"""
del r_coord
phi_coord = phi_coord.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_phi_coord', phi_coord)
"""
del phi_coord
theta_coord = theta_coord.reshape((len(z),1))
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_theta_coord', theta_coord)
"""
del theta_coord
bass_decam = bass_decam.reshape((len(z),1))  #True decam, False bass
""" NOTE changed here
np.save('/cluster/scratch/ntheiler/ucat/ucat_sorted_app_mag_20220526_pin500mpc20488_6_2_zlt05_123_m_bass_decam', bass_decam)
"""
del bass_decam
host_sub_index = host_sub_index.reshape((len(z),1))

del z
del halo_mass
del x_coord
del y_coord
del z_coord
del host_sub_index

print('output prepared')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))

print('output saved')
print("--- %s seconds ---" % round(time.clock() - start_time, 2))


# ---------------------------------------
# PRINTOUTS

print(' ')
print('pinocchio_halo_subhalo_catalog: ' + str(pinocchio_halo_subhalo_catalog))
#print('outfile_galaxies: ' + str(outfile_galaxies))
print('output: abs_mag, z, app_mag_g, app_mag_r, app_mag_z, blue_red (blue=1), r_coord, phi_coord, theta_coord, bass_decam (decam=1)')
print(' ')
print('Number of galaxies in total: ' + str(n_gal))
print('Number of red galaxies: ' + str(n_gal - np.count_nonzero(blue_red)))
print('Number of red galaxies with non-nan magnitudes: ' + str(np.count_nonzero(~np.isnan(abs_mag) & (blue_red < 0.5))))
print('Number of blue galaxies: ' + str(np.count_nonzero(blue_red)))
print('Number of blue galaxies with non-nan magnitudes: ' + str(np.count_nonzero(~np.isnan(abs_mag) & (blue_red > 0.5))))
