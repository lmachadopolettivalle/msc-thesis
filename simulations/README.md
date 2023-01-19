# Simulations

## Step 1: PINOCCHIO

To run PINOCCHIO, execute `./submit_pinocchio.sh`. This will change into the appropriate directory (where we wish to store the output data, i.e. in `scratch`), submit the job using Slurm, and change back into the original directory. As a result, the input files can live at `home/` while the outputs are saved in `scratch/`.

Reminders:
- Modify `submit_pinocchio.sh` to `cd` into the appropriate output directory
- Modify `parameter_file` with specific particle count, box size, redshift range, and lightcone aperture angle
- Modify `plc_desired_redshifts` to contain the desired redshifts at which output files will be saved (Necessary when modigying redshift range in `parameter_file`)
- Modify `sbatch_pinocchio_script.sh` for appropriate run time, memory requested, and job name


## Step 2: Extract subhalos from PINOCCHIO outputs

Execute `./submit_subhalos.sh` to submit the subhalo code job. This will use the generated PINOCCHIO halo catalog and merger histories to extract subhalos and save them into a new series of output files. The generated subhalos can then be used for the SHAM model.

Prerequisites:
- Completed PINOCCHIO run

Reminders:
- Modify `sbatch_pinocchio_script.sh` for appropriate run time, memory requested, and job name
- Modify `compute_subhalos_from_pinocchio.py` for appropriate:
	- Particle count
	- Footprint mask input file
	- Directory containing PINOCCHIO outputs
	- PINOCCHIO Slurm output file
	- Run flag used in PINOCCHIO

## Step 3: Sample luminosity functions for red and blue galaxies (Can be run independently from Step 2)

`NOTE: this only needs to happen once! No need to rerun it for every single SHAM model.`

Execute `./submit_sample_only.sh` to submit the sample job. This will sample galaxies from a luminosity function, generating a catalog of many red and blue galaxies (separately) with their redshifts and absolute magnitudes.

Prerequisites:
- Completed PINOCCHIO run (purely because it needs the cosmology used, which is read in from the PINOCCHIO output file)

## Steps 4 + 5: Construct 2D histogram for SHAM + Perform SHAM assignment

The next two steps will often happen one after the other. As a result, there is a script designed to submit both Steps as Slurm jobs.

This is especially useful when exploring many sets of SHAM parameters, as this script reduces the risk of losing track of which parameters were used for any given run.

To submit both steps together, execute `submit_2D_sham.py` with command line arguments corresponding to the parameters desired for the run.

Or, to submit many sets of parameters at once, execute `submit_many_runs.sh` instead.

Prerequisites:
- Completed PINOCCHIO run
- Completed subhalo extraction code
- Completed sampling of galaxies from luminosity function

Reminders:
- If using `submit_many_runs.sh`, modify `submit_many_runs.sh` with the correct sets of desired SHAM parameters to vary. Remember to set default values for the other parameters as well.

### Submitting Steps 4 and 5 individually
If, instead of submitting Steps 4 and 5 together with `submit_many_runs.sh`, you desire to submit either step individually, follow these steps.

### Step 4: Construct 2D histogram for SHAM

Execute `./submit_2D_hist.sh` to submit the 2D histogram job. It reads the extracted subhalos from PINOCCHIO, and divides them into redshift and mass bins, which are then used in the SHAM assignment algorithm.

Prerequisites:
- Completed PINOCCHIO run
- Completed subhalo extraction code

### Step 5: Perform the SHAM assignment

Execute `./submit_sham_interpolation.sh` to submit the SHAM code job. This performs the following tasks:

- Compute apparent magnitudes for the sampled galaxies, using instrument-specific filters (e.g. BASS/MzLS and DECam)
- Assign galaxies to halos and subhalos with SHAM, using the 2D histograms from a previous step

Prerequisites:
- Completed PINOCCHIO run
- Completed subhalo extraction code
- Completed sampling of galaxies from luminosity function
- Completed generation of 2D histograms

## Evaluation: magnitude and color histograms, and 2PCF
