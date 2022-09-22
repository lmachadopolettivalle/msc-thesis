from desitarget.cuts import select_targets
import numpy as np

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

FILES = ["sweep-040p025-050p030.fits", "sweep-260p025-270p030.fits", "sweep-010p000-020p005.fits", "sweep-020p025-030p030.fits", "sweep-240p025-250p030.fits"]

# Select BGS targets from sweep files
targets, infiles = select_targets(
    FILES,
    numproc=1,
    tcnames=["BGS"],
    survey="main",
    backup=False,
    return_infiles=True
)

# Load RA, DEC from BGS targets
RA_index = -1
DEC_index = -1

# Find indices for desired columns
for i, column in enumerate(targets.dtype.descr):
    if column[0] == "RA":
        RA_index = i
    if column[0] == "DEC":
        DEC_index = i

    if (RA_index != -1) and (DEC_index != -1):
        break

RAs = [i[RA_index] for i in targets]
DECs = [i[DEC_index] for i in targets]
N = len(targets)

# Get random RAs and DECs
# TODO Need to use the DR9 randoms instead of np.random
rand_N = 2000
random_RAs = 360 * np.random.rand(rand_N)
random_DECs = 180 * np.random.rand(rand_N) - 90

# Compute correlation function
nbins = 10
bins = np.linspace(0.1, 10.0, nbins + 1)
nthreads = 1

DD_counts = DDtheta_mocks(1, nthreads, bins, RAs, DECs)
RR_counts = DDtheta_mocks(1, nthreads, bins, random_RAs, random_DECs)
DR_counts = DDtheta_mocks(0, nthreads, bins, RAs, DECs, RA2=random_RAs, DEC2=random_DECs)

wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
DD_counts, DR_counts,
DR_counts, RR_counts)

print(wtheta)
