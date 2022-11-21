#!/bin/bash
#SBATCH --mem-per-cpu=4000
#SBATCH -n 1032
#SBATCH -N 9-32
#SBATCH -t 01:00:00
#SBATCH -J pinocchio_job_512
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_512
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/pinocchio_error_512
#SBATCH --mail-type=FAIL

mpirun -np 1032 /cluster/home/lmachado/PINOCCHIO/Pinocchio/src/pinocchio.x /cluster/home/lmachado/msc-thesis/simulations/parameter_file
