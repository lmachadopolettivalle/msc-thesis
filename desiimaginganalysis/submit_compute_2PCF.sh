#!/bin/bash
#SBATCH --mem-per-cpu=6400
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 04:00:00
#SBATCH -J luis_2PCF_job
#SBATCH -o luis_2PCF_output
#SBATCH -e luis_2PCF_error
#SBATCH --mail-type=ALL

/cluster/home/lmachado/venv/bin/python -u compute_2PCF.py
