#!/bin/bash
#SBATCH --mem-per-cpu=3900
#SBATCH -n 100
#SBATCH -N 1
#SBATCH -t 00:45:00
#SBATCH --mail-type=FAIL

/cluster/home/lmachado/venv/bin/python -u construct_2D_hists.py --run_id ${run_id} --particle_count_pinocchio ${particle_count_pinocchio} --region ${region}
