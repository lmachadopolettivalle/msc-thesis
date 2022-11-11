# Simulations

## PINOCCHIO

To run PINOCCHIO, execute `./submit_pinocchio.sh`. This will change into the appropriate directory (where we wish to store the output data, i.e. in `scratch`), submit the job using Slurm, and change back into the original directory. As a result, the input files can live at `home/` while the outputs are saved in `scratch/`.

To visualize the PINOCCHIO outputs, execute `python PlotExample.py`. This will read the outputs from `scratch` and show the Mass Function and the positions of the generated halos.
- Make sure to change `PlotExample.py` to use the same run name as the one used in the `parameter_file`.
