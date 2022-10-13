# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np
from desitarget.targets import bgs_mask

# Parameters for magnitude computation
CLIP_FLUX = 1e-16

REGION = "north"

BANDS = ["g", "r", "z"]

# Parameters for plotting
bin_count = 90

LINEWIDTH = 2

# Color choices
blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"

plt.rcParams['font.size'] = '12'
bright_plot_color = blue
faint_plot_color = red

# If True, add magnitudes without extinction correction
# on top of magnitude histograms
# If False, only plot extinction corrected magnitudes
compare_extinction_correction = False

# Load target data
with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_BGS_TARGET.npy", "rb") as f:
    bgs_targets = np.load(f)

fluxes = {}
mw_transmissions = {}
for band in BANDS:
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_FLUX_{band.upper()}.npy", "rb") as f:
        fluxes[band] = np.load(f)
    with open(f"/cluster/scratch/lmachado/DataProducts/targets/{REGION}/targets_MW_TRANSMISSION_{band.upper()}.npy", "rb") as g:
        mw_transmissions[band] = np.load(g)

# Separate BGS Bright and BGS Faint targets
bright_targets = np.array([i for i, bgs_target in enumerate(bgs_targets) if "BGS_BRIGHT" in bgs_mask.names(bgs_target)])
faint_targets = np.array([i for i, bgs_target in enumerate(bgs_targets) if "BGS_FAINT" in bgs_mask.names(bgs_target)])

def mag_from_flux(targets, band="r", extinction_correction=True):
    # Extinction-corrected magnitude in g, r, z bands

    flux = fluxes[band][targets]
    if extinction_correction:
        mw_transmission = mw_transmissions[band][targets]
        mags = 22.5 - 2.5 * np.log10((flux / mw_transmission).clip(CLIP_FLUX))
    else:
        mags = 22.5 - 2.5 * np.log10(flux.clip(CLIP_FLUX))

    return mags

def remove_spurious_mag_objects(mags_dict):
    # Modify mags_dict to remove entries from all its mags,
    # if mag in at least one of them is due to clipping of the flux.
    # Explanation: due to the clip(1e-16), there will be a few objects (<10 or so)
    # with magnitudes == 62.5.
    # These objects can be ignored.

    clip_mag = 22.5 - 2.5 * np.log10(CLIP_FLUX)

    r_idx = mags_dict["r"] < clip_mag
    g_idx = mags_dict["g"] < clip_mag
    z_idx = mags_dict["z"] < clip_mag

    all_idx = r_idx & g_idx & z_idx

    clipped_mags_dict = {
        k: v[all_idx]
        for k, v in mags_dict.items()
    }

    return clipped_mags_dict


# Dicts to hold magnitudes in each band
bright_mags = {}
faint_mags = {}

if compare_extinction_correction:
    bright_mags_not_extinction_corrected = {}
    faint_mags_not_extinction_corrected = {}

# Plot histograms of magnitudes in 3 bands
for band in BANDS:
    bright_mags[band] = np.array(mag_from_flux(bright_targets, band=band, extinction_correction=True))
    faint_mags[band] = np.array(mag_from_flux(faint_targets, band=band, extinction_correction=True))

# Remove unhelpful objects, which have spurious flux values
bright_mags = remove_spurious_mag_objects(bright_mags)
faint_mags = remove_spurious_mag_objects(faint_mags)

print("Number of bright objects:", len(bright_mags["r"]))
print("Number of faint objects:", len(faint_mags["r"]))

# Determine binning for each color band
bins = {
    band: np.linspace(
        min(bright_mags[band]),
        max(faint_mags[band]),
        bin_count,
    )
    for band in BANDS
}

# Plot histograms
for band in BANDS:
    if compare_extinction_correction:
        bright_label = "BGS Bright (E.C.)"
        faint_label = "BGS Faint (E.C.)"
    else:
        bright_label = "BGS Bright"
        faint_label = "BGS Faint"

    plt.hist(bright_mags[band], bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="step", label=bright_label, color=bright_plot_color)
    plt.hist(faint_mags[band], bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="step", label=faint_label, color=faint_plot_color)

    if compare_extinction_correction:
        bright_mags_not_extinction_corrected[band] = np.array(mag_from_flux(bright_targets, band=band, extinction_correction=False))
        faint_mags_not_extinction_corrected[band] = np.array(mag_from_flux(faint_targets, band=band, extinction_correction=False))

        plt.hist(bright_mags_not_extinction_corrected[band], bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="stepfilled", alpha=0.5, label="BGS Bright (Not E.C.)", color=bright_plot_color)
        plt.hist(faint_mags_not_extinction_corrected[band], bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="stepfilled", alpha=0.5, label="BGS Faint (Not E.C.)", color=faint_plot_color)

    plt.xlim([14, 22])
    #plt.ylim([0, 2])
    if compare_extinction_correction:
        plt.title(f"{band} mag, with and without extinction correction (E.C.)")
        plt.xlabel(f"{band} magnitude")
    else:
        plt.title(f"{band} mag distribution")
        plt.xlabel(f"{band} magnitude, extinction-corrected")

    plt.ylabel("PDF")
    plt.legend(loc="upper left")
    plt.grid()
    if compare_extinction_correction:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{band}mag_hist_compare_extinction.pdf")
    else:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{band}mag_hist.pdf")

    #plt.show()
    plt.clf()

# Plot histograms of colors
bright_gminusr = bright_mags["g"] - bright_mags["r"]
faint_gminusr = faint_mags["g"] - faint_mags["r"]

bright_rminusz = bright_mags["r"] - bright_mags["z"]
faint_rminusz = faint_mags["r"] - faint_mags["z"]

# G - R Color Histogram
for density in {True, False}:
    plt.hist(bright_gminusr, bins=bin_count, linewidth=LINEWIDTH, density=density, histtype="step", label="BGS Bright", color=bright_plot_color)
    plt.hist(faint_gminusr, bins=bin_count, linewidth=LINEWIDTH, density=density, histtype="step", label="BGS Faint", color=faint_plot_color)

    plt.xlim([-0.5, 2.5])
    plt.xlabel("g - r color")
    if density:
        plt.ylabel("PDF")
    else:
        plt.ylabel("Count")
    plt.legend(loc="upper left")
    plt.grid()
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/gminusr_{'density' if density else 'nodensity'}_hist.pdf")
    #plt.show()
    plt.clf()

    # R - Z Color Histogram
    plt.hist(bright_rminusz, bins=bin_count, linewidth=LINEWIDTH, density=density, histtype="step", label="BGS Bright", color=bright_plot_color)
    plt.hist(faint_rminusz, bins=bin_count, linewidth=LINEWIDTH, density=density, histtype="step", label="BGS Faint", color=faint_plot_color)

    plt.xlim([-0.5, 2.5])
    plt.xlabel("r - z color")
    if density:
        plt.ylabel("PDF")
    else:
        plt.ylabel("Count")
    plt.legend(loc="upper left")
    plt.grid()
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/rminusz_{'density' if density else 'nodensity'}_hist.pdf")
    #plt.show()
    plt.clf()

#####
# G-R vs. R-Z color plot
# Since this is a scatter plot,
# cannot handle plotting all targets.
# Instead, sample a subset of the targets to be plotted.
sample_number = 10000
bright_sample_idx = np.random.choice(np.arange(len(bright_rminusz)), sample_number, replace=False)
faint_sample_idx = np.random.choice(np.arange(len(faint_rminusz)), sample_number, replace=False)

plt.scatter(bright_rminusz[bright_sample_idx], bright_gminusr[bright_sample_idx], label="BGS Bright", s=4, alpha=0.5, c=bright_plot_color)
plt.scatter(faint_rminusz[faint_sample_idx], faint_gminusr[faint_sample_idx], label="BGS Faint", s=4, alpha=0.5, c=faint_plot_color)

plt.xlabel("r - z color")
plt.ylabel("g - r color")
plt.xlim([-1.5, 3.5])
plt.ylim([-1.5, 3.5])
plt.legend(loc="upper right")
plt.grid()
plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/color.pdf")
#plt.show()
plt.clf()
