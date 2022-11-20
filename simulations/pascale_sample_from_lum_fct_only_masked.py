# Drawing Galaxies from Luminosity Functions

# UCat only, just get redshifts and absolute magnitudes, and red vs. blue
# No apparent magnitudes, no matching to PINOCCHIO
# Pixarea calculated from a mask

# author: Pascale Berner
# first written: 09.11.2022
# last adapted: 16.11.2022
# partially copied from: sample_from_lum_fct_only.py

# ----------------------------------------------------
# IMPOTRS
# -----------------------------------------------------
from ucat import galaxy_sampling_util # need ucatenv2 or test
from ucat import io_util
import h5py
import scipy.integrate
import PyCosmo
import numpy as np

# ----------------------------------------------------
# PARAMETERS TO BE SET
# -----------------------------------------------------
z_max = 0.75
m_max = -12

pixarea = 0 # set pixarea to the desired value, if no mask is provided, otherwise leave it at 0
infile_footprint = '/cluster/work/refregier/bernerp/DES/DES_Y1/y1a1_gold_1.0.2_wide_footprint_4096.fit'
NSIDE = 4096
dec_min = ???
dec_max = ???

outfile_dir = './'
outfile_z = 'sampled_z_desy1_20221116_1.npy'
outfile_absmag = 'sampled_absmag_desy1_20221116_1.npy'
outfile_redblue = 'sampled_redblue_desy1_20221116_1.npy'

omega_b=0.045
omega_m=0.3
omega_l_in=0.7
pk_norm=0.81
n_param=0.961

# ----------------------------------------------------
# OTHER PARAMETERS
# -----------------------------------------------------
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

m_star_par_blue = (lum_fct_m_star_blue_slope, 
                   lum_fct_m_star_blue_intcpt)
m_star_par_red = (lum_fct_m_star_red_slope,
                  lum_fct_m_star_red_intcpt)
phi_star_par_blue = (lum_fct_phi_star_blue_amp,
                     lum_fct_phi_star_blue_exp)
phi_star_par_red = (lum_fct_phi_star_red_amp,
                    lum_fct_phi_star_red_exp)

# ----------------------------------------------------
# SET PYCOSMO COSMOLOGY
# -----------------------------------------------------
cosmo = PyCosmo.Cosmo()
cosmo.set(omega_b=omega_b, omega_m=omega_m, omega_l_in=omega_l_in, pk_norm=pk_norm, n=n_param)
cosmo.print_params()

# ----------------------------------------------------
# LOAD MASK, CALCULATE PIXAREA
# -----------------------------------------------------
if pixarea == 0:
	m_footprint = hp.fitsfunc.read_map(infile_footprint)
	m_dec_range = ????
	n_pix = np.sum((m_footprint > 0) & (m_dec_range > 0))
	pixarea = n_pix / hp.nside2npix(NSIDE) * 4. * np.pi

# ----------------------------------------------------
# SAMPLE GALAXIES
# -----------------------------------------------------
n_gal_cal_blue = galaxy_sampling_util.NumGalCalculator(z_max, m_max, lum_fct_alpha_blue, m_star_par_blue, phi_star_par_blue, cosmo, pixarea)
n_gal_cal_red = galaxy_sampling_util.NumGalCalculator(z_max, m_max, lum_fct_alpha_red, m_star_par_red, phi_star_par_red, cosmo, pixarea)

n_gal_blue = n_gal_cal_blue()
n_gal_red = n_gal_cal_red()

z_mabs_sampler_blue = galaxy_sampling_util.RedshiftAbsMagSampler(lum_fct_z_res, z_max, lum_fct_m_res, m_max, lum_fct_alpha_blue,
																m_star_par_blue, phi_star_par_blue, cosmo)

z_mabs_sampler_red = galaxy_sampling_util.RedshiftAbsMagSampler(lum_fct_z_res, z_max, lum_fct_m_res, m_max, lum_fct_alpha_red, 
																m_star_par_red, phi_star_par_red, cosmo)

z_blue, abs_mag_blue = z_mabs_sampler_blue(n_gal_blue, m_max, lum_fct_alpha_blue, m_star_par_blue)
z_red, abs_mag_red = z_mabs_sampler_red(n_gal_red, m_max, lum_fct_alpha_red, m_star_par_red)
index_blue = np.ones(len(z_blue)) # 1 = blue
index_red = np.zeros(len(z_red)) # 0 = red

n_tot = n_gal_red + n_gal_blue
z_tot = np.concatenate((z_red, z_blue))
abs_mag_tot = np.concatenate((abs_mag_red, abs_mag_blue))
index_tot = np.concatenate((index_red, index_blue))

# ----------------------------------------------------
# SAVE GALAXIES TO NPY FILES
# -----------------------------------------------------
np.save(outfile_dir+outfile_z, z_tot)
np.save(outfile_dir+outfile_absmag, absmag_tot)
np.save(outfile_dir+outfile_redblue, index_tot)

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
print('This is a pure UCat run, until absolute magnitudes and redshifs are drawn. No apparent magnitudes etc. are calculated, no matching to PINOCCHIO is done.') 


