#!/bin/bash
#SBATCH --mem-per-cpu=3900
#SBATCH -n 128
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -J subhalos_job_fullsky_2048_fullsky
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/subhalos_output_fullsky_2048_fullsky
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/subhalos_error_fullsky_2048_fullsky
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u compute_subhalos_from_pinocchio.py
