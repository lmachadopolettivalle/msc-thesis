#!/bin/bash
#SBATCH --mem-per-cpu=40000
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J subhalos_job
#SBATCH -o /cluster/home/lmachado/msc-thesis/simulations/subhalos_output
#SBATCH -e /cluster/home/lmachado/msc-thesis/simulations/subhalos_error
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python compute_subhalos_pascale.py
