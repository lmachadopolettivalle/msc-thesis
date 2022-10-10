# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np
from desitarget.targets import bgs_mask

blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"

plt.rcParams['font.size'] = '12'
bright_plot_color = blue
faint_plot_color = yellow

bin_count = 75

# If True, add magnitudes without extinction correction
# on top of magnitude histograms
# If False, only plot extinction corrected magnitudes
compare_extinction_correction = False

# Load target data
with open("/cluster/scratch/lmachado/DataProducts/targets/targets_BGS_TARGET.npy", "rb") as f:
    bgs_targets = np.load(f)

fluxes = {}
mw_transmissions = {}
for band in {"g", "r", "z"}:
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/targets_FLUX_{band.upper()}.npy", "rb") as f:
        fluxes[band] = np.load(f)
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/targets_MW_TRANSMISSION_{band.upper()}.npy", "rb") as g:
        mw_transmissions[band] = np.load(g)

# Separate BGS Bright and BGS Faint targets
bright_targets = np.array([i for i, bgs_target in enumerate(bgs_targets) if "BGS_BRIGHT" in bgs_mask.names(bgs_target)])
faint_targets = np.array([i for i, bgs_target in enumerate(bgs_targets) if "BGS_FAINT" in bgs_mask.names(bgs_target)])

def mag_from_flux(targets, band="r", extinction_correction=True):
    # Extinction-corrected magnitude in g, r, z bands
    flux = fluxes[band][targets]
    if extinction_correction:
        mw_transmission = mw_transmissions[band][targets]
        return 22.5 - 2.5 * np.log10(flux / mw_transmission)
    else:
        return 22.5 - 2.5 * np.log10(flux)

# Dicts to hold magnitudes in each band
bright_mags = {}
faint_mags = {}

if compare_extinction_correction:
    bright_mags_not_extinction_corrected = {}
    faint_mags_not_extinction_corrected = {}

# Plot histograms of magnitudes in 3 bands
for band in ["g", "r", "z"]:
    bright_mags[band] = np.array(mag_from_flux(bright_targets, band=band, extinction_correction=True))
    faint_mags[band] = np.array(mag_from_flux(faint_targets, band=band, extinction_correction=True))

    if compare_extinction_correction:
        bright_label = "BGS Bright (E.C.)"
        faint_label = "BGS Faint (E.C.)"
    else:
        bright_label = "BGS Bright"
        faint_label = "BGS Faint"

    plt.hist(bright_mags[band], bins=bin_count, density=True, histtype="step", label=bright_label, color=bright_plot_color)
    plt.hist(faint_mags[band], bins=bin_count, density=True, histtype="step", label=faint_label, color=faint_plot_color)

    if compare_extinction_correction:
        bright_mags_not_extinction_corrected[band] = np.array(mag_from_flux(bright_targets, band=band, extinction_correction=False))
        faint_mags_not_extinction_corrected[band] = np.array(mag_from_flux(faint_targets, band=band, extinction_correction=False))

        plt.hist(bright_mags_not_extinction_corrected[band], bins=bin_count, density=True, histtype="stepfilled", alpha=0.5, label="BGS Bright (Not E.C.)", color=bright_plot_color)
        plt.hist(faint_mags_not_extinction_corrected[band], bins=bin_count, density=True, histtype="stepfilled", alpha=0.5, label="BGS Faint (Not E.C.)", color=faint_plot_color)

    plt.xlim([14, 22])
    plt.ylim([0, 2])
    if compare_extinction_correction:
        plt.title(f"{band} mag, with and without extinction correction (E.C.)")
        plt.xlabel(f"{band} magnitude")
    else:
        plt.title(f"{band} mag distribution")
        plt.xlabel(f"{band} magnitude, extinction-corrected")

    plt.ylabel("PDF")
    plt.legend()
    plt.grid()
    if compare_extinction_correction:
        plt.savefig(f"{band}mag_hist_compare_extinction.pdf")
    else:
        plt.savefig(f"{band}mag_hist.pdf")

    plt.show()

# Plot histograms of colors
bright_gminusr = bright_mags["g"] - bright_mags["r"]
faint_gminusr = faint_mags["g"] - faint_mags["r"]

bright_rminusz = bright_mags["r"] - bright_mags["z"]
faint_rminusz = faint_mags["r"] - faint_mags["z"]

# G - R Color Histogram
plt.hist(bright_gminusr, bins=bin_count, density=True, histtype="step", label="BGS Bright", color=bright_plot_color)
plt.hist(faint_gminusr, bins=bin_count, density=True, histtype="step", label="BGS Faint", color=faint_plot_color)

plt.xlim([-0.5, 2.5])
plt.xlabel("g - r color")
plt.ylabel("PDF")
plt.legend()
plt.grid()
plt.savefig("gminusr_hist.pdf")
plt.show()

# R - Z Color Histogram
plt.hist(bright_rminusz, bins=bin_count, density=True, histtype="step", label="BGS Bright", color=bright_plot_color)
plt.hist(faint_rminusz, bins=bin_count, density=True, histtype="step", label="BGS Faint", color=faint_plot_color)

plt.xlim([-0.5, 2.5])
plt.xlabel("r - z color")
plt.ylabel("PDF")
plt.legend()
plt.grid()
plt.savefig("rminusz_hist.pdf")
plt.show()

# G-R vs. R-Z color plot
plt.scatter(bright_rminusz, bright_gminusr, label="BGS Bright", s=2, c=bright_plot_color)
plt.scatter(faint_rminusz, faint_gminusr, label="BGS Faint", s=2, c=faint_plot_color)

plt.xlabel("r - z color")
plt.ylabel("g - r color")
plt.legend()
plt.grid()
plt.savefig("color.pdf")
plt.show()
