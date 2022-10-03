# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np
from desitarget.targets import bgs_mask

BINS = 100

# Load target data
with open("targets.npy", "rb") as f:
    data = np.load(f)

# Separate BGS Bright and BGS Faint targets
bright_data = [i for i in data if "BGS_BRIGHT" in bgs_mask.names(i["BGS_TARGET"])]
faint_data = [i for i in data if "BGS_FAINT" in bgs_mask.names(i["BGS_TARGET"])]

def mag_from_flux(target, band="r"):
    # Extinction-corrected magnitude in g, r, z bands
    band = band.upper()
    flux = target[f"FLUX_{band}"]
    mw_transmission = target[f"MW_TRANSMISSION_{band}"]
    return 22.5 - 2.5 * np.log10(flux / mw_transmission)

# Dicts to hold magnitudes in each band
bright_mags = {}
faint_mags = {}

# Plot histograms of magnitudes in 3 bands
for band in ["g", "r", "z"]:
    bright_mags[band] = np.array([mag_from_flux(i, band=band) for i in bright_data])
    faint_mags[band] = np.array([mag_from_flux(i, band=band) for i in faint_data])

    plt.hist(bright_mags[band], bins=BINS, histtype="step", label="BGS Bright")
    plt.hist(faint_mags[band], bins=BINS, histtype="step", label="BGS Faint")

    plt.xlabel(f"{band} magnitude, extinction-corrected")
    plt.ylabel("Count")
    plt.legend()
    plt.grid()
    plt.savefig(f"{band}mag_hist.pdf")
    plt.show()

# Plot histograms of colors
bright_gminusr = bright_mags["g"] - bright_mags["r"]
faint_gminusr = faint_mags["g"] - faint_mags["r"]

bright_rminusz = bright_mags["r"] - bright_mags["z"]
faint_rminusz = faint_mags["r"] - faint_mags["z"]

# G - R Color Histogram
plt.hist(bright_gminusr, bins=BINS, histtype="step", label="BGS Bright")
plt.hist(faint_gminusr, bins=BINS, histtype="step", label="BGS Faint")

plt.xlabel("g - r color")
plt.ylabel("Count")
plt.legend()
plt.grid()
plt.savefig("gminusr_hist.pdf")
plt.show()

# R - Z Color Histogram
plt.hist(bright_rminusz, bins=BINS, histtype="step", label="BGS Bright")
plt.hist(faint_rminusz, bins=BINS, histtype="step", label="BGS Faint")

plt.xlabel("r - z color")
plt.ylabel("Count")
plt.legend()
plt.grid()
plt.savefig("rminusz_hist.pdf")
plt.show()


# G-R vs. R-Z color plot
plt.scatter(bright_rminusz, bright_gminusr, label="BGS Bright", s=2)
plt.scatter(faint_rminusz, faint_gminusr, label="BGS Faint", s=2)

plt.xlabel("r - z color")
plt.ylabel("g - r color")
plt.legend()
plt.grid()
plt.savefig("color.pdf")
plt.show()
