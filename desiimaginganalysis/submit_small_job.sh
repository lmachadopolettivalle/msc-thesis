#!/bin/bash
#SBATCH --mem-per-cpu=2000
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH -J luis_job
#SBATCH -o luis_output
#SBATCH -e luis_error

/cluster/home/lmachado/venv/bin/python $1
