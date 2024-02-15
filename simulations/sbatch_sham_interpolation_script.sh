#!/bin/bash
#SBATCH --mem-per-cpu=4096
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL

/cluster/home/lmachado/venv/bin/python -u sham_interpolation.py --run_id ${run_id} --particle_count_pinocchio ${particle_count_pinocchio} --region ${region} --seed ${seed}
