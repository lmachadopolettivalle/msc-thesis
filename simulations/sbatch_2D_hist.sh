#!/bin/bash
#SBATCH --mem-per-cpu=4096
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 00:20:00
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u construct_2D_hists.py --run_id ${run_id} --particle_count_pinocchio ${particle_count_pinocchio}
