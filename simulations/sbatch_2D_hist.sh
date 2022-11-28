#!/bin/bash
#SBATCH --mem-per-cpu=4096
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 00:20:00
#SBATCH -J 2D_hist_job_2048
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/2D_hist_output_2048
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/2D_hist_error_2048
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u construct_2D_hists.py --num_z_bins 150 --num_mass_bins 30
