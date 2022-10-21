import numpy as np

NUMBER_OF_POINTS = 30000000

RNG = np.random.default_rng(42)

U, V = np.array(list(zip(*RNG.random((NUMBER_OF_POINTS, 2)))))

RA = np.degrees(2 * np.pi * U)
DEC = np.degrees(np.pi / 2 - np.arccos(2 * V - 1))


with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_RA.npy", "wb") as f:
    np.save(f, RA)
with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_DEC.npy", "wb") as f:
    np.save(f, DEC)
