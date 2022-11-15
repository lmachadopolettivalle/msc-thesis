bsub -n 10 -R "span[ptile=2]" -R "rusage[mem=3500]" mpirun -np 10 /cluster/home/lmachado/PINOCCHIO/Pinocchio/src/pinocchio.x /cluster/home/lmachado/msc-thesis/simulations/parameter_file
