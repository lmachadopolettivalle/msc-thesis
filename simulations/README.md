# Simulations

## Step 1: PINOCCHIO

To run PINOCCHIO, execute `./submit_pinocchio.sh`. This will change into the appropriate directory (where we wish to store the output data, i.e. in `scratch`), submit the job using Slurm, and change back into the original directory. As a result, the input files can live at `home/` while the outputs are saved in `scratch/`.

To visualize the PINOCCHIO outputs, execute `python PlotExample.py`. This will read the outputs from `scratch` and show the Mass Function and the positions of the generated halos. NOTE: make sure to modify this script to read in multiple files, if you requested PINOCCHIO to store its output in many files, per the `NumFiles` setting.

- Make sure to change `PlotExample.py` to use the same run name as the one used in the `parameter_file`.

## Step 2: Extract subhalos from PINOCCHIO outputs

Execute `./submit_subhalos.sh` to submit the subhalo code job. This will use the generated PINOCCHIO halo catalog and merger histories to extract subhalos and save them into a new series of output files. The generated subhalos can then be used for the SHAM model.

Prerequisites:
- Completed PINOCCHIO run

## Step 3: Construct 2D histogram for SHAM

Execute `./submit_2D_hist.sh` to submit the 2D histogram job. It reads the extracted subhalos from PINOCCHIO, and divides them into redshift and mass bins, which are then used in the SHAM assignment algorithm.

Prerequisites:
- Completed PINOCCHIO run
- Completed subhalo extraction code

## Step 4: Sample luminosity functions for red and blue galaxies (Can be run independently from Step 2)

Execute `./submit_sample_only.sh` to submit the sample job. This will sample galaxies from a luminosity function, generating a catalog of many red and blue galaxies (separately) with their redshifts and absolute magnitudes.

Prerequisites:
- Completed PINOCCHIO run (purely because it needs the cosmology used, which is read in from the PINOCCHIO output file)

## Step 5: Perform the SHAM assignment

Execute `./submit_sham_interpolation.sh` to submit the SHAM code job. This performs the following tasks:

- Compute apparent magnitudes for the sampled galaxies, using instrument-specific filters (e.g. BASS/MzLS and DECam)
- Assign galaxies to halos and subhalos with SHAM, using the 2D histograms from a previous step

Prerequisites:
- Completed PINOCCHIO run
- Completed subhalo extraction code
- Completed sampling of galaxies from luminosity function
- Completed generation of 2D histograms

## Evaluation: magnitude and color histograms, and 2PCF
