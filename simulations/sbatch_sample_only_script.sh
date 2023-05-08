#!/bin/bash
#SBATCH --mem-per-cpu=800
#SBATCH -n 128
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -J sample_only_job_fullsky_2048
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/sample_only_output_fullsky_2048
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/sample_only_error_fullsky_2048
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u sample_from_lumfct_masked.py
