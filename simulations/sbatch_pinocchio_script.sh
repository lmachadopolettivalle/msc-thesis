#!/bin/bash
#SBATCH --mem-per-cpu=4000
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J pinocchio_job
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/pinocchio_output
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/pinocchio_error
#SBATCH --mail-type=FAIL

mpirun -np 4 /cluster/home/lmachado/PINOCCHIO/Pinocchio/src/pinocchio.x /cluster/home/lmachado/msc-thesis/simulations/parameter_file
