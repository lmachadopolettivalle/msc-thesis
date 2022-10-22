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

## Compute HEALPIx pixel IDs for selected targets

Run `compute_target_healpix.py` to use the targets' RA and DEC values to determine the HEALPix pixel IDs where the targets lie on the sky.

Prerequisites:
- After running `save_filtered_targets_to_files.py`
- Modify `compute_target_healpix.py` to use `REGION = "north"` or `REGION = "south"`, as well as to choose the desired NSIDE

## Count number of generated targets (optional)

If interested, can compute number of objects in the original sweep files to compare against the resulting number of selected BGS targets. For this, use `count_targets.py`.

Prerequisites:
- Run target selection
- Modify `count_targets.py` to use `REGION = "north"` or `REGION = "south"`, as well as to choose the desired NSIDE

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
- Modify `mask.py` to use `REGION = "north"` or `REGION = "south"`, as well as to choose the desired NSIDE

## Visualization: Plot 1-point functions (histograms)

Run `plot_histograms.py` (using `Slurm`, e.g. `sbatch submit_small_job.sh plot_histograms.py`) to generate magnitude and color histograms based on the selected BGS targets.

Prerequisites:
- Run target selection
- Modify `plot_histograms.py` to use `REGION = "north"` or `REGION = "south"`, as well as to choose the desired NSIDE

Results:
- Save several PDF files (under `images/`) with the histogram plots

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

Then, use `plot_2PCF.py` to display the computed corelation function.

Prerequisites:
- Run target selection
- Run code to draw random points on the sphere
