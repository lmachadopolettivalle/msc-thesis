from matplotlib import pyplot as plt
import numpy as np

REGION = "north"

with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_bins.npy", "rb") as f:
    bins = np.load(f)
with open(f"/cluster/scratch/lmachado/DataProducts/2PCF/{REGION}_2PCF_wtheta.npy", "rb") as f:
    wtheta = np.load(f)

plt.plot(bins, wtheta)

plt.xscale("log")
plt.yscale("log")

plt.show()
