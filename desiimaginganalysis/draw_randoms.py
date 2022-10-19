import numpy as np

NUMBER_OF_POINTS = 100

RNG = np.random.default_rng()

U, V = np.array(list(zip(*RNG.random((NUMBER_OF_POINTS, 2)))))

RA = np.degrees(2 * np.pi * U)
DEC = np.degrees(np.pi / 2 - np.arccos(2 * V - 1))


