#!/bin/bash
#SBATCH --mem-per-cpu=4096
#SBATCH -n 128
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J sham_interpolation_job_2048
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/sham_interpolation_output_2048
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/sham_interpolation_error_2048
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u sham_interpolation.py
