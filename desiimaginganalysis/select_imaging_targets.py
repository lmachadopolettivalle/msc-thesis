from desitarget.cuts import select_targets
from desitarget.targets import desi_mask, bgs_mask
import os
import time

REGION = "south"

PATH_TO_SWEEP_FILES = f"/cluster/scratch/lmachado/DESIImaging/dr9/{REGION}/sweeps/"

FILENAMES = os.listdir(PATH_TO_SWEEP_FILES)

FILES = [f"{PATH_TO_SWEEP_FILES}/{filename}" for filename in FILENAMES]

start = time.time()
print("Starting at", start)

# Select BGS targets from sweep files
targets, infiles = select_targets(
    FILES,
    numproc=16,
    tcnames=["BGS"],
    survey="main",
    backup=False,
    return_infiles=True,
    resolvetargs=True,
)

end = time.time()
print("Finished selection at", end)
print("Time elapsed in seconds:", end-start)

# Load RA, DEC from BGS targets
RAs = [i["RA"] for i in targets]
DECs = [i["DEC"] for i in targets]
N = len(targets)
