# Description of the program:
# extract mergers from merger history of PINOCCHIO,
# load the halo catalog (for multiple files of the halo lightcone),
# apply a healpix mask, calculate the surviving subhalos, and
# save the resulting halo-subhalo catalog (in separate files)

# Author: Pascale Berner
# Co-Author: Luis Machado
# first written: 04.11.2022
# last adapted: 24.11.2022
# partially copied from: halo_subhalo_from_plc_hist_scatterinpos_new_6.py

# ------------------------------------------
# IMPORTS
# ------------------------------------------


print("Importing required libraries...")

from collections import defaultdict
from halotools.empirical_models import NFWProfile # ATTENTION: needs hdf5 and python/3.6.0!!!
import healpy as hp
import numpy as np
import os
import pandas as pd
import re
from tqdm import tqdm

np.random.seed(42)

print("Done importing libraries.")

# Initialize NFWProfile model
nfw_model = NFWProfile()

# -----------------------------------------------------
# SPECIFICATIONS FOR SIMULATION RUNS
# -----------------------------------------------------
# THIS IS A SECTION THAT NEEDS TO BE ADAPTED EACH TIME!
# -----------------------------------------------------

# Number of files to be processed in this run
# Ideally this could be the total number of files.
# However, due to the large memory and time requirements,
# we need to split this run into many batches,
# which means we only process a few files at a time.
NUMBER_OF_FILES_TO_BE_PROCESSED = 24

# Cube root of number of particles used.
# This is present in the paths to different input files used in this script.
particle_count_pinocchio = 2048

# Directory containing PINOCCHIO outputs.
# This is where the subhalo catalog will also be saved at the end of this script.
dirname = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/"

pinocchio_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_{particle_count_pinocchio}" # Path to SLURM output from PINOCCHIO, which contains many useful details on the run

RUN_FLAG = "luis" # Corresponds to "RunFlag" in PINOCCHIO parameter file. Name of the run being analyzed.

cosmology_file = f'pinocchio.{RUN_FLAG}.cosmology.out'
plc_file = f'pinocchio.{RUN_FLAG}.plc.out'
history_file = f'pinocchio.{RUN_FLAG}.histories.out'

# Create output directory
outfile_dir = f"/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/luis_runs/{particle_count_pinocchio}cubed/halo_subhalo_plc/"
if os.path.isdir(outfile_dir):
    print(f"{outfile_dir} directory already exists.")
else:
    print(f"Creating new output directory, {outfile_dir} ...")
    os.mkdir(outfile_dir)
    print("Created output directory successfully.")

outfile_halos = 'pinocchio_masked_halos_subhalos_plc'

# ------------------------------------------
# MASKING CONFIGURATION
# ------------------------------------------

# Decide which HEALPix ordering to use
# Set NEST = False for RING, NEST = True for NEST
NEST = True

# Provide mask used to filter regions in which to find subhalos
# Mask should be a HEALPix map,
# with desired pixel IDs having a value > MASK_CUTOFF_VALUE,
# and undesired pixel IDs having a value <= MASK_CUTOFF_VALUE
# E.g. an array with 0s and 1s could use a cutoff value of 0
MASK_CUTOFF_VALUE = 0

# We support an additional masking based on a declination range.
# If this feature is not desired, set dec_min = -90 and dec_max = 90
dec_min = -90
dec_max = 90

# Path to file containing mask array as a HEALPix map
# If no mask is desired, set the filename to None, and the full-sky data will be included in the analysis
infile_footprint = "/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy"

# Set a default NSIDE and NPIX and a full-sky mask, in case there is no valid input mask
NSIDE = 1
NPIX = hp.pixelfunc.nside2npix(NSIDE)
m_footprint = np.ones(NPIX)

# Load pixel mask, if any is provided
if infile_footprint is not None:
    print(f"Loading mask from {infile_footprint}...")

    if infile_footprint.endswith(".fit"):
        m_footprint = hp.fitsfunc.read_map(infile_footprint, nest=NEST)
    elif infile_footprint.endswith(".npy"):
        with open(infile_footprint, "rb") as f:
            m_footprint = np.load(f)
    else:
        print(f"WARNING: Invalid file format for HEALPix mask file {infile_footprint}. Using a full-sky mask instead.")

    # Determine NSIDE, NPIX from provided map
    NSIDE = hp.get_nside(m_footprint)
    NPIX = len(m_footprint)
    assert hp.pixelfunc.nside2npix(NSIDE) == NPIX

    print("Finished loading mask.")
else:
    print("No mask file provided. Will use default full-sky mask instead.")

# ------------------------------------------
# BEGINNING OF SUBHALO CODE.
# NO NEED TO MODIFY BELOW THIS LINE.
# ------------------------------------------

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
    m_part = float(
        re.search("Particle Mass \(Msun/h\)\s+(.+)\n", text).groups()[0]
    )
    min_num_part = int(
        re.search("MinHaloMass \(particles\)\s+(\d+)\n", text).groups()[0]
    )
    num_rep = int(
        re.search("The box will be replicated (\d+) times to construct the PLC\n", text).groups()[0]
    )
    z_max, z_min = (
        float(i)
        for i in re.search("The Past Light Cone will be reconstruct from z=(.+) to z=(.+)\n", text).groups()
    )
    omega_l = float(
        re.search("OmegaLambda\s+(.+)\n", text).groups()[0]
    )
    omega_m = float(
        re.search("Omega0\s+(.+)\n", text).groups()[0]
    )
    H0 = 100 * float(
        re.search("Hubble100\s+(.+)\n", text).groups()[0]
    )

    return (
        num_files,
        m_part,
        min_num_part,
        num_rep,
        z_min,
        z_max,
        omega_l,
        omega_m,
        H0,
    )

num_files, m_part, min_num_part, num_rep, z_min, z_max, omega_l, omega_m, H0 = read_pinocchio_config_details()

# Determine file endings based on number of catalog files
if num_files == 1:
    file_ending = [""]
else:
    file_ending = [f".{i}" for i in range(num_files)]

# print directory etc.
print('dirname = ' + dirname)
print('cosmology_file = ' + cosmology_file)
print('plf_file = ' + plc_file)
print('history_file = ' + history_file)
print('outfile_halos = ' + outfile_halos)
print('outfile_dir = ' + outfile_dir)
print('m_part = ' + str(m_part))
print('min_num_part = ' + str(min_num_part))
print('num_rep = ' + str(num_rep))
print('num_files = ' + str(num_files))
print('z_min = ' + str(z_min))
print('z_max = ' + str(z_max))
print('NSIDE = ' + str(NSIDE))
print('dec_min = ' + str(dec_min))
print('dec_max = ' + str(dec_max))
print('infile_footprint = ' + str(infile_footprint))

# -----------------------------------------------------
# SPECIFICATIONS FOR THE USED COSMOLOGY
# -----------------------------------------------------
H0_1 = 1./(H0/(3.09e19))/(3.15576e16) # inverse of hubble constant in Gyr
t_dyn_0 = 0.1*H0_1 # dynamical time, using H0 instead of H(z)

# -----------------------------------------------------
# Preparation for sampling for eta (using pdf eqn 12 of Simon Birrers paper 1401.3162)
# -----------------------------------------------------

eta_x = np.linspace(0,1,101)
pdf = eta_x**(1.2)*(1.-eta_x)**(1.2) # calculating pdf values
cdf = np.cumsum(pdf) # calculating cdf values
cdf = cdf/cdf[-1] # normalize cdf
cdf_min = cdf[eta_x == 0.2]
#cdf_max = cdf[eta_x == 0.8] # used earlier
cdf_max = cdf[eta_x == 1.0] # does not matter much

print('definitions done')

# -----------------------------------------------------
# DEFINITIONS OF FUNCTIONS
# -----------------------------------------------------

# calculate time difference for two redshifts
# eqn. (3.99) in Mo, Bosch and White
def delta_time_from_z_z2array(z1, z2): # z2 being an array, z1 a float
	t1 = H0_1*2./(3*np.sqrt(omega_l))*np.log((np.sqrt(omega_l*(1+z1)**(-3.))+np.sqrt(omega_l*(1+z1)**(-3.)+omega_m))/np.sqrt(omega_m))
	t1_array = t1*np.ones(len(z2))
	return t1_array-H0_1*2./(3*np.sqrt(omega_l))*np.log((np.sqrt(omega_l*(np.ones(len(z2))+z2)**(-3.))+np.sqrt(omega_l*(np.ones(len(z2))+z2)**(-3.)+omega_m))/np.sqrt(omega_m))

def t_dyn(z):
	return t_dyn_0/np.sqrt(omega_l + omega_m*(1+z)**3)

# preparation for position of subhalo within the host
# sampling using an Isothermal profile, according to arXiv:0402160
def cum_num(rad):
    #rad: normalised radius, rad = r/(0.37*r_vir)
    return rad - np.arctan(rad)
rad_x = np.linspace(0., 1./0.37, 1000) # reference array for sampling the normalised radius rad
cum_N = cum_num(rad_x)
Nmax = np.max(cum_N)
Nmin = 0.

# loop over haloid, get subs out of dict2
# input of loop: l=number of halos
# output of loop: l=number of halos and subhalos combined
#@profile
def loop_dict(l):
    for i in range(l):
        subs = np.array(subs_dict[haloid[i]])
        if len(subs) != 0:
            factor = 1.0
            b_val = 0.92*growth_rate[i]
            tdf = (0.9*t_dyn_0*0.216*np.ones(len(subs)))*((m_ratio[subs])**b_val * factor)/np.log(np.ones(len(subs))+m_ratio[subs])*np.exp(1.9*eta[subs])
            delta_t = delta_time_from_z_z2array(redshift[i], acc_red[subs]) # delta_t can be negative like this
            subhalos = subs[(delta_t < tdf) & (delta_t >= 0) & (m_ratio_temp[subs] < 100.)]
            n = len(subhalos)
            if n > 0:
                delta_t_subhalos = delta_t[(delta_t < tdf) & (delta_t >= 0) & (m_ratio_temp[subs] < 100.)] # for saving
                tdf_subhalos = tdf[(delta_t < tdf) & (delta_t >= 0) & (m_ratio_temp[subs] < 100.)] # for saving
                N_sample = np.random.uniform(Nmin, Nmax, n)
                x_sample = np.array([rad_x[cum_N == min(cum_N[(cum_N - yy) > 0])][0] for yy in N_sample])
                rad_i = x_sample * 0.37 * nfw_model.halo_mass_to_halo_radius(M[i]) # assumption: the given radius is approximately the virial radius
                x_i = np.reshape(x_rand_prep[l:l+n]*rad_i+X[i], (n,1))
                y_i = np.reshape(y_rand_prep[l:l+n]*rad_i+Y[i], (n,1))
                z_i = np.reshape(z_rand_prep[l:l+n]*rad_i+Z[i], (n,1))

                halo_subhalo_array[l:l+n] = np.concatenate((np.reshape(groupid[subhalos], (n,1)), np.reshape(num_part[subhalos], (n,1))*m_part, redshift[i]*np.ones((n,1)), x_i, y_i, z_i, np.zeros((n,1)), M[i]*np.ones((n,1)), np.zeros((n,1)), np.reshape(delta_t_subhalos, (n,1)), np.reshape(tdf_subhalos, (n,1)), haloid[i]*np.ones((n,1)) ), axis=1)
                halo_subhalo_array[i,6] = M[i] - m_part*np.sum(num_part[subhalos])
                l = l+n
    return l

# -----------------------------------------------------
# LOAD SCALE FACTOR AND GROWTH RATE FROM COSMOLOGY FILE
# -----------------------------------------------------

# Get linear growth rate for each redshift
data = pd.read_csv(dirname+cosmology_file, sep='\s+', lineterminator='\n', header=None, index_col=None, skipinitialspace=True).values
scale_factor_sorted = data[:, 0]
growth_rate_sorted = data[:, 2]
# preparation for interpolation later
redshift_sorted = 1./scale_factor_sorted - 1.

# -----------------------------------------------------
# IMPORT MERGER HISTORY
# -----------------------------------------------------
print("Loading merger history...")

groupid = np.array([])
treeind = np.array([])
merged_with = np.array([])
num_part = np.array([])
num_part_merged_with = np.array([])
acc_red = np.array([])
for i in range(num_files):
    data = pd.read_csv(dirname+history_file+file_ending[i], sep='\s+', lineterminator='\n', header=None, index_col=None, skiprows=16, comment='#').values
    groupid_temp = data[:, 0]
    treeind_temp = data[:, 1]
    merged_with_temp = data[:, 3]
    num_part_temp = data[:, 4]
    num_part_merged_with_temp = data[:, 5]
    acc_red_temp = data[:, 6]

    groupid = np.concatenate((groupid, groupid_temp), axis=None)
    treeind = np.concatenate((treeind, treeind_temp), axis=None)
    merged_with = np.concatenate((merged_with, merged_with_temp), axis=None)
    num_part = np.concatenate((num_part, num_part_temp), axis=None)
    num_part_merged_with = np.concatenate((num_part_merged_with, num_part_merged_with_temp), axis=None)
    acc_red = np.concatenate((acc_red, acc_red_temp), axis=None)

n_groups = len(groupid)

# -----------------------------------------------------
# CREATE HISTORY DICTIONARY
# -----------------------------------------------------

# using history, create dictionary containing haloid in 1st column and position in history of main host at z=0 in 2nd column
# -> dict1
mainhost_groupid_dict = defaultdict(list)

# using history and dict1, create dictionary containing haloid in 1st column and list of groupids of groups that merged into this halo (all subgroups) in 2nd column
# -> dict2
subs_dict = defaultdict(list)

mainhost_id = -1
num_group_members = -1

# recursion to get halo one rank up and write current halo as a subhalo of that into dict2
def write_host_to_dict(sub_id_ind, host_linking):
	if host_linking != -1:
		one_up_id = int(mainhost_id + (host_linking % num_group_members))
		subs_dict[int(groupid[one_up_id])].append(int(sub_id_ind))
		write_host_to_dict(sub_id_ind, merged_with[one_up_id])

# loop once through the history to prepare dict1 and dict2
for id_ind, id_val in enumerate(groupid):
	host_link = merged_with[id_ind] # id within group of halo higher up
	if host_link == -1: # check if it is a host
		mainhost_id = id_ind # adapt the current host
		num_group_members = treeind[id_ind]
	else:
		write_host_to_dict(id_ind, host_link)
	mainhost_groupid_dict[id_val].append(mainhost_id) # save id in history of current host

# -----------------------------------------------------
# PREPARATIONS OUTSIDE THE LOOP
# -----------------------------------------------------

# mass ratio between halos and subhalos, capped at 40
m_ratio = 1.*num_part_merged_with/num_part # M_host/M_satellite at accretion
m_ratio_temp = m_ratio
m_ratio = np.minimum(m_ratio, 40.*np.ones(len(m_ratio)))

cdf_sample = np.random.uniform(cdf_min, cdf_max, n_groups)
eta = np.array([eta_x[cdf == min(cdf[(cdf - yy) > 0])][0] for yy in cdf_sample])
# rcrv: R_c/R_vir, parameter for the orbital energy
rcrv = np.random.uniform(0.1, 1., n_groups)

print("Finished loading history.")

# -----------------------------------------------------
# PREPARE FOR POSITIONS OF SURVIVING SUBHALOS WITHIN HOSTS
# -----------------------------------------------------

# sampling in direction
rand_ra = np.random.uniform(-np.pi, np.pi, n_groups*num_rep)
rand_sindec = np.random.uniform(-1., 1., n_groups*num_rep)
rand_dec = np.arcsin(rand_sindec)
x_rand_prep = np.sin(np.pi/2+rand_dec)*np.cos(rand_ra)
y_rand_prep = np.sin(np.pi/2+rand_dec)*np.sin(rand_ra)
z_rand_prep = np.cos(np.pi/2+rand_dec)
#ignoring that massive subhalos should be distributed differently than light ones...

# -----------------------------------------------------
# LOOP OVER LIGHTCONE FILES
# -----------------------------------------------------

# counter over all halos: n_halos
# counter over all halos and subhalos: n_all
# counter for halos and subhalos within one file: l
n_halos = 0
n_subhalos = 0
n_all = 0

# Keep track of how many files have been processed.
# Stop when we reach the limit set by the user above,
# or when we have completed all files.
count_files_processed = 0

for i in range(num_files):
    # ------------------------------
    # CHECK IF WE HAVE PROCESSED ENOUGH FILES,
    # AND STOP IF WE HAVE.
    # ------------------------------
    if count_files_processed >= NUMBER_OF_FILES_TO_BE_PROCESSED:
        print(f"Finished processing the minimum requested number of files: {NUMBER_OF_FILES_TO_BE_PROCESSED}")
        print(f"Actually processed {count_files_processed} files.")
        print("Stopping this job...")
        break

    # If we have not processed the requested number of files,
    # begin next file
    print(f"Starting computation for file ending {i}...")
    # ------------------------------
    # CHECK IF FILE HAS ALREADY BEEN PROCESSED.
    # IF SO, CONTINUE TO NEXT FILE
    # ------------------------------
    input_halo_filename = dirname + plc_file + file_ending[i]
    output_subhalo_filename = dirname + outfile_halos + file_ending[i] + ".txt"

    if os.path.exists(output_subhalo_filename):
        print(f"Files have already been created for file ending {i}, skipping this index...")
        continue

    # ------------------------------
    # LOAD ONE LIGHTCONE FILE
    # ------------------------------
    data = pd.read_csv(input_halo_filename, sep='\s+', lineterminator='\n', header=None, index_col=None, skipinitialspace=True).values
    haloid = data[:, 0]
    redshift = data[:, 1]
    X = data[:, 2]
    Y = data[:, 3]
    Z = data[:, 4]
    M = data[:, 8]
    # ------------------------------
    # APPLY PIXEL MASK and DEC LIMITS
    # ------------------------------
    # calculate RA and DEC
    r_calc_pin = np.sqrt(X**2 + Y**2 + Z**2) # radial distance
    theta_calc_pin = np.arccos(Z/r_calc_pin) # theta
    phi_calc_pin = np.arctan2(Y, X) # phi

    # Note that phi is in range [-pi, pi], but for healpy, must be in range [0, 360 degrees]
    phi_calc_pin[phi_calc_pin < 0] += 2 * np.pi

    ra = np.degrees(phi_calc_pin) # RA
    dec = np.degrees(np.pi/2 - theta_calc_pin)

    # APPLY MASK AND DEC LIMITS
    pix_halos = hp.ang2pix(NSIDE, ra, dec, lonlat=True, nest=NEST)
    spatial_mask = (m_footprint[pix_halos] > MASK_CUTOFF_VALUE) & (dec >= dec_min) & (dec <= dec_max)

    M = M[spatial_mask]
    haloid = haloid[spatial_mask]
    redshift = redshift[spatial_mask]
    X = X[spatial_mask]
    Y = Y[spatial_mask]
    Z = Z[spatial_mask]
    # ------------------------------
    # PREPARE HALO-SUBHALO ARRAY
    # ------------------------------
    l = len(haloid)
    halo_subhalo_array = np.zeros((l*3,12))
    halo_subhalo_array[:l,0] = haloid # column 0: ID of halo itself
    halo_subhalo_array[:l,1] = M # column 1: Mass of halo or subhalo itself
    halo_subhalo_array[:l,2] = redshift # column 2: Redshift of halo (for subhalos: redshift of host)
    halo_subhalo_array[:l,3] = X # column 3, 4, 5: position of halo
    halo_subhalo_array[:l,4] = Y
    halo_subhalo_array[:l,5] = Z
    halo_subhalo_array[:l,6] = M # column 6: Mass of halo minus mass of its subhalos
    # column 7: Mass of the halos host (for hosts: 0)
    halo_subhalo_array[:l,8] = np.ones(l) # column 8: 1 for hosts, 0 for subhalos
    # column 9: delta_t (for hosts: 0)
    # column 10: tdf (for hosts: 0)
    # column 11: ID of host (for hosts: 0)
    n_halos += l
    # ------------------------------
    # CALCULATE GROWTH RATE FROM Z
    # ------------------------------
    growth_rate = np.interp(redshift, redshift_sorted, growth_rate_sorted)
    # ------------------------------
    # RUN INNER LOOP OVER HALOS
    # ------------------------------
    l = loop_dict(l)
    n_subhalos += (l - len(haloid))
    n_all += l
    # ------------------------------
    # SAVE HALO-SUBHALO ARRAY TO FILE
    # ------------------------------
    # order by descending mass, along the 1st column, only first l rows
    halo_subhalo_array_sorted = halo_subhalo_array[halo_subhalo_array[:,1].argsort()][::-1][:l,:]

    with open(output_subhalo_filename, 'w') as f:
        np.savetxt(f, halo_subhalo_array_sorted, fmt='%d %.8e %.8f %.8e %.8e %.8e %.8e %.8e %d %.8e %.8e %d')

    count_files_processed += 1
    print(f"Finished processing files for file ending {i}.")

# -----------------------------------------------------
# SOME PRINT OUTPUT
# -----------------------------------------------------
print('n_halos = ' + str(n_halos))
print('n_subhalos = ' + str(n_subhalos))
print('n_all = ' + str(n_all))
print('outfile_halos: groupid, group mass [Msun], redshift of host, x, y, z of host [Mpc/h], mass of host minus masses of subhalos, mass of host, host?, delta_t, tdf, halo id of host')

print("Job complete.")
