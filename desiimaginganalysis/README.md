# DESI Legacy Survey Analysis Pipeline

## Download data

 We use the `sweep` files available at https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/ and https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/

Run `download_sweep_data.py` at the IPA cluster. Then, `rsync` the downloaded sweepfiles to Euler, where we perform all the computations.

Also, can modify `download_sweep_data.py` to download `footprint` files instead of `sweep` files. This is useful for the creation of a footprint HEALPix mask, which is optional.

Prerequisites: None.

Result:
- Download sweepfiles to the cluster.

## Perform BGS Target Selection

Run `save_filtered_targets_to_files.py` in Euler (using `Slurm`, e.g. `sbatch submit_job.sh save_filtered_targets_to_files.py`) to read the sweepfiles and obtain selected BGS target objects.

Prerequisites:
- Download sweepfiles to Euler
- Modify `select_imaging_targets.py` to use `REGION = "north"` or `REGION = "south"`

Result:
- Create several `target_*.npy` files, one for each data column.

## Useful counts of targets and mask areas (optional)

If interested, can compute a few useful counts for the selected targets and mask areas, as well as the number of objects in the original sweep files to compare against the resulting number of selected BGS targets. For this, use `useful_target_counts.py`.

Prerequisites:
- Run target selection
- Modify `useful_target_counts.py` to use `REGION = "north"` or `REGION = "south"`, as well as to choose the desired NSIDE

## Create mask based on footprint area (optional)

Run `create_footprint_mask.py` to create a mask based on the intended Legacy Survey bricks that were planned to be observed.

Note that this is a superset of the sky areas that actually contain targets. During the comparison against simulations, as well as during the computation of the correlation function, we do NOT use this mask, but instead use the mask based on the actual selected target positions.

Prerequisites:
- Download footprint files to the cluster

Results:
- Generate footprint mask files.

## Create mask based on selected target positions

Run `mask.py` to generate HEALPix masks for both north and south regions. This generated mask lies entirely within the target area, guaranteeing no issues with boundary pixels while also minimizing area loss.

Prerequisites:
- Run target selection for both north and south regions
- Modify `mask.py` to choose the desired NSIDE

Results:
- Generate HEALPix mask files for BASS/MzLS, DECaLS-NGC, and DECaLS-SGC based on target data
- Create mollview plot of targets and masks

# Load target data (helper file)

The file `load_processed_target_data.py` offers functionality to load the selected BGS targets and return their key values in a Python dictionary.

If requested, this functionality will also apply the mask for the desired subset of regions (BASS/MzLS, DECaLS-NGC, DECaLS-SGC) to the target data, only returning targets within the chosen regions.

Prerequisites:
- Running target selection
- Running script to create masks

## Visualization: Plot 1-point functions (histograms)

Run `plot_histograms.py` (using `Slurm`, e.g. `sbatch submit_small_job.sh plot_histograms.py`) to generate magnitude and color histograms based on the selected BGS targets.

Prerequisites:
- Run target selection
- Modify `plot_histograms.py` to choose the desired NSIDE

Results:
- Save several PDF files (under `images/`) with the histogram plots

## Visualization: color scatter plot separated by target's MORPHTYPE

To see color distributions of each MORPHTYPE selected, use `plot_morphtypes_colors.py` to create a scatter plot with colors of targets separated by MORPHTYPE.

Prerequisites:
- Run target selection

Results:
- Save a PDF file (under `images/`) with the color scatter plot

## Visualization: plot targets on the sky, as well as footprints

Use `plot_sky_map.py` to display the selected targets as well as the footprint on a Mollview plot.

Prerequisites:
- Run target selection
- Run HEALPix code to compute pixel IDs of targets
- Run code to generate footprint based on Legacy Survey files

Results:
- Save PDF file (under `images/`) with Mollview plot of the target and footprint areas

## Draw randoms for correlation function

Use `draw_randoms.py` (with `Slurm`) to generate large sets of RAs and DECs and save them to a file.

Prerequisites: None.

Results:
- Save `npy` files with random RA and DEC values.

## Compute angular correlation function

Use `compute_2PCF.py` (with `Slurm`) to compute the angular correlation function of the selected targets, using the LS estimator.

Then, use `plot_2PCF.py` to display the computed correlation function.

Prerequisites:
- Run target selection
- Run code to draw random points on the sphere
- Modify code to choose which REGION (BASS/MzLS, DECaLS-NGC, DECaLS-SGC) to use when computing 2PCF
