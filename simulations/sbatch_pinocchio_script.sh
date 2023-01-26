#!/bin/bash
#SBATCH --constraint=ib
#SBATCH --mem-per-cpu=6000
#SBATCH -n 1032
#SBATCH -N 9-32
#SBATCH -t 01:00:00
#SBATCH -J pinocchio_job_deep_2048
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_deep_2048
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/pinocchio_error_deep_2048
#SBATCH --mail-type=FAIL

mpirun -np 1032 /cluster/home/lmachado/PINOCCHIO/Pinocchio/src/pinocchio.x /cluster/home/lmachado/msc-thesis/simulations/parameter_file
