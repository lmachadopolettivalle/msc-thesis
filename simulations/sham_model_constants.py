import numpy as np
# -----------------------------------------------------
# UCAT / UFIG PARAMETERS
# -----------------------------------------------------

template_coeff_z1_blue = 1                                     # Redshift z1>0 for blue galaxies
template_coeff_z1_red = 1                                      # Redshift z1>0 for red galaxies
magnitude_calculation = 'table'                                # The way magnitudes are calculated
lum_fct_filter_band = 'GenericBessel_B'                        # Filter band in which the luminosity function is valid

maps_remote_dir = 'ufig_res/maps/'
filters_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/filters_collection.h5'
templates_int_tables_file_name = '/cluster/work/refregier/bernerp/DESI_BGS_filters/ufig_res/sed_integrals__template_spectra_BlantonRoweis07.h5'


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

