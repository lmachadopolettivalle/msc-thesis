#!/bin/bash
#SBATCH --mem-per-cpu=30G
#SBATCH -n 16
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH -J luis_big_job
#SBATCH -o luis_big_output
#SBATCH -e luis_big_error
#SBATCH --mail-type=END,FAIL

/cluster/home/lmachado/venv/bin/python $1
