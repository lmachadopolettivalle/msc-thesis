#!/bin/bash
#SBATCH --mem-per-cpu=6400
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J compute_2PCF_job
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/compute_2PCF_output
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/compute_2PCF_error
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u compute_2PCF.py --region "fullsky" --run_id 162
# 140, 157, 158, 159, 160, 161, 162
