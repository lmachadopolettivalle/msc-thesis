# Given a run configuration and a particle count,
# submit the sbatch jobs to create a 2D histogram,
# and to run the SHAM algorithm.

import argparse
import re
import subprocess

from manage_parameter_space import get_id_of_run, store_new_run

# Obtain relevant command-line arguments
parser = argparse.ArgumentParser()

parser.add_argument("--region", type=str, required=True)
parser.add_argument("--num_z_bins", type=int, required=True)
parser.add_argument("--num_mass_bins", type=int, required=True)
parser.add_argument("--mass_cut", type=float, required=True)
parser.add_argument("--quenching_time", type=float, required=True)
parser.add_argument("--particle_count_pinocchio", type=int, required=True)

args = parser.parse_args()

region = args.region
num_z_bins = args.num_z_bins
num_mass_bins = args.num_mass_bins
mass_cut = args.mass_cut
quenching_time = args.quenching_time
particle_count_pinocchio = args.particle_count_pinocchio

# Save this run to the database, and then obtain its stored ID
store_new_run(num_z_bins=num_z_bins, num_mass_bins=num_mass_bins, mass_cut=mass_cut, quenching_time=quenching_time)

run_id = get_id_of_run(num_z_bins=num_z_bins, num_mass_bins=num_mass_bins, mass_cut=mass_cut, quenching_time=quenching_time)

# Submit sbatch jobs using this run ID
submit_message_2D = subprocess.run(
    [
        "sbatch",
        f"--export=ALL,run_id={run_id},particle_count_pinocchio={particle_count_pinocchio}",
        f"-J2D_hist_job_deep_{particle_count_pinocchio}_{run_id}",
        f"-o/cluster/home/lmachado/msc-thesis/simulations/2D_hist_deep_{particle_count_pinocchio}_{run_id}_output",
        f"-e/cluster/home/lmachado/msc-thesis/simulations/2D_hist_deep_{particle_count_pinocchio}_{run_id}_error",
        "sbatch_2D_hist.sh",
    ],
    capture_output=True,
)

stdout_message_2D = submit_message_2D.stdout.decode("utf-8")
job_id_2D = int(
    re.search("Submitted batch job (\d+)", stdout_message_2D).groups()[0]
)

submit_message_sham = subprocess.run(
    [
        "sbatch",
        f"--dependency=afterok:{job_id_2D}",
        f"--export=ALL,run_id={run_id},particle_count_pinocchio={particle_count_pinocchio},region={region}",
        f"-Jsham_int_job_deep_{particle_count_pinocchio}_{run_id}",
        f"-o/cluster/home/lmachado/msc-thesis/simulations/sham_int_job_deep_{particle_count_pinocchio}_{run_id}_output",
        f"-e/cluster/home/lmachado/msc-thesis/simulations/sham_int_job_deep_{particle_count_pinocchio}_{run_id}_error",
        "sbatch_sham_interpolation_script.sh",
    ],
)
