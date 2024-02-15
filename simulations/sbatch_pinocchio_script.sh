#!/bin/bash
#SBATCH --constraint=ib
#SBATCH --mem-per-cpu=6000
#SBATCH -n 1032
#SBATCH -N 9-32
#SBATCH -t 02:00:00
#SBATCH -J pinocchio_job_fullsky_2048_666666
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/pinocchio_output_fullsky_2048_666666
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/pinocchio_error_fullsky_2048_666666
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lmachado@phys.ethz.ch

mpirun -np 1032 /cluster/home/lmachado/PINOCCHIO/Pinocchio/src/pinocchio.x /cluster/home/lmachado/msc-thesis/simulations/parameter_file
