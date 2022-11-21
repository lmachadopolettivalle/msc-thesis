#!/bin/bash
#SBATCH --mem-per-cpu=40000
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J sample_only_job_512
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/sample_only_output_512
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/sample_only_error_512
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python sample_from_lumfct_masked.py
