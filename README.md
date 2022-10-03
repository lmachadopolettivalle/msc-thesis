# MSC Thesis

## Legacy Imaging Survey, DR9

- Data files: sweep files
- Pipeline:
	- Download files from the web: `download_sweep_data.py`
	- Select targets from sweep files using `desitarget`: `select_imaging_targets.py`
	- Save selected targets in `.npy` format: `save_filtered_targets_to_files.py`
	- Compute 2PCF: `compute_2pcf.py`

- Useful code:
	- Generate mask/footprint from survey brick files: `generate_mask_footprint.py`
	- Test memory cuts from `desitarget`: `test_memory_cut_from_target_selection.py`
