from select_imaging_targets import *

from matplotlib import pyplot as plt
import numpy as np

import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

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

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\theta$")
plt.ylabel(r"w($\theta$)")
plt.scatter(
    bins[:-1],
    wtheta
)
plt.savefig("test2PCF.pdf")
