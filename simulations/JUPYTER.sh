#!/bin/bash
#SBATCH --mem-per-cpu=2000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH --output JUPYTER.out


# Steps:
# 1. Submit this job with sbatch
# 2. Once approved, see the output file. It will have a ssh command. Run that ssh in a new, local terminal. That will not print any logs, but it will be working!
# 3. Then, check again for the output file. It will contain a link to 127.0.0.1/... - copy this link into your browser.

let ipnport=4200
echo ipnport=$ipnport
ipnip=$(hostname -i)

echo "ssh -N -L $ipnport:$ipnip:$ipnport lmachado@euler.ethz.ch"

export XDG_RUNTIME_DIR=""

source ~/venv/bin/activate

jupyter notebook --no-browser --ip 0.0.0.0 --port=$ipnport
