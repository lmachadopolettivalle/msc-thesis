#!/bin/bash
#SBATCH --mem-per-cpu=40000
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J sham_interpolation_job
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/sham_interpolation_output
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/sham_interpolation_error
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python sham_interpolation.py
