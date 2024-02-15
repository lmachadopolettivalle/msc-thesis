PINOCCHIO_BASE = "/cluster/scratch/lmachado/PINOCCHIO_OUTPUTS/"

HISTOGRAMS = "2D_histograms"
INTERPOLATION_OUTPUTS = "interpolation_outputs"
TWOPCF = "2PCF"
OUTPUTS_SAMPLING = "outputs_sampling"

BASS_MzLS = "BASS-MzLS"
DECaLS_NGC = "DECaLS-NGC"
DECaLS_SGC = "DECaLS-SGC"
FULLSKY = "fullsky"

MASK_FILES = {
    BASS_MzLS: "/cluster/scratch/lmachado/DataProducts/masks/BASS_MzLS_mask.npy",
    DECaLS_NGC: "/cluster/scratch/lmachado/DataProducts/masks/DECaLS_NGC_mask.npy",
    DECaLS_SGC: "/cluster/scratch/lmachado/DataProducts/masks/DECaLS_SGC_mask.npy",
    FULLSKY: "/cluster/scratch/lmachado/DataProducts/masks/fullsky_mask.npy",
}

# NOTE: most PINOCCHIO runs will be stored under the seedXXXXXX path. But the original PINOCCHIO run was not. This is why we have the if...else block within this function
def pinocchio_halo_files_path(particle_count=2048, z_depth=0.5, pinocchio_region="north", seed=None):
    if seed is None:
        SEED_PATH = ""
    else:
        SEED_PATH = f"seed{seed}"

    return f"{PINOCCHIO_BASE}/{SEED_PATH}/{particle_count}cubed/z{z_depth:.1f}/{pinocchio_region}/"

def pinocchio_subhalo_files_path(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", seed=None):
    halo_path = pinocchio_halo_files_path(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, seed=seed)

    return f"{halo_path}/{DESI_region}/"

def outputs_sampling_path(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", seed=None):
    subhalo_path = pinocchio_subhalo_files_path(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, DESI_region=DESI_region, seed=seed)

    return f"{subhalo_path}/{OUTPUTS_SAMPLING}/"

def path_run(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", run_id=100, seed=None):
    subhalo_path = pinocchio_subhalo_files_path(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, DESI_region=DESI_region, seed=seed)

    return f"{subhalo_path}/{run_id}/"

def path_2PCF(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", run_id=100, seed=None):
    run = path_run(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, DESI_region=DESI_region, run_id=run_id, seed=seed)

    return f"{run}/{TWOPCF}/"

def path_histograms(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", run_id=100, seed=None):
    run = path_run(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, DESI_region=DESI_region, run_id=run_id, seed=seed)

    return f"{run}/{HISTOGRAMS}/"

def path_interpolation(particle_count=2048, z_depth=0.5, pinocchio_region="north", DESI_region="BASS-MzLS", run_id=100, seed=None):
    run = path_run(particle_count=particle_count, z_depth=z_depth, pinocchio_region=pinocchio_region, DESI_region=DESI_region, run_id=run_id, seed=seed)

    return f"{run}/{INTERPOLATION_OUTPUTS}/"
