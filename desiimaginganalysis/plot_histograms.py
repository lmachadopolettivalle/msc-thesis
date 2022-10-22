# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np

from load_processed_target_data import load_processed_target_data, BANDS, BGS_BRIGHT, BGS_FAINT

REGION = "north"

# If True, use magnitudes with extinction correction
# If False, show magnitudes without extinction correction
extinction_correction = True

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

# Load targets
targets = load_processed_target_data(region=REGION, extinction_correction=extinction_correction)

bright_targets = np.where(targets["BGS_TARGET"] == BGS_BRIGHT)[0]
faint_targets = np.where(targets["BGS_TARGET"] == BGS_FAINT)[0]

# Plot histograms of magnitudes in 3 bands

print("Number of bright objects:", len(targets["MAG_R"][bright_targets]))
print("Number of faint objects:", len(targets["MAG_R"][faint_targets]))

# Determine binning for each color band
bins = {
    band: np.linspace(
        min(targets[f"MAG_{band.upper()}"][bright_targets]),
        max(targets[f"MAG_{band.upper()}"][faint_targets]),
        bin_count,
    )
    for band in BANDS
}

# Plot histograms
for band in BANDS:
    if extinction_correction:
        bright_label = "BGS Bright"
        faint_label = "BGS Faint"
    else:
        bright_label = "BGS Bright (E.C.)"
        faint_label = "BGS Faint (E.C.)"

    plt.hist(
        targets[f"MAG_{band.upper()}"][bright_targets],
        bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="step", label=bright_label, color=bright_plot_color
    )
    plt.hist(
        targets[f"MAG_{band.upper()}"][faint_targets],
        bins=bins[band], linewidth=LINEWIDTH, density=False, histtype="step", label=faint_label, color=faint_plot_color
    )

    plt.xlim([14, 22])
    #plt.ylim([0, 2])
    if extinction_correction:
        plt.title(f"{band} mag distribution")
        plt.xlabel(f"{band} magnitude, extinction-corrected")
    else:
        plt.title(f"{band} mag, with and without extinction correction (E.C.)")
        plt.xlabel(f"{band} magnitude")

    plt.ylabel("Count")
    plt.legend(loc="upper left")
    plt.grid()
    if extinction_correction:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_{band}mag_hist.pdf")
    else:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_{band}mag_hist_no_extinction.pdf")

    #plt.show()
    plt.clf()

# Plot histograms of colors
bright_gminusr = targets["MAG_G"][bright_targets] - targets["MAG_R"][bright_targets]
faint_gminusr = targets["MAG_G"][faint_targets] - targets["MAG_R"][faint_targets]

bright_rminusz = targets["MAG_R"][bright_targets] - targets["MAG_Z"][bright_targets]
faint_rminusz = targets["MAG_R"][faint_targets] - targets["MAG_Z"][faint_targets]

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
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_gminusr_{'density' if density else 'nodensity'}_hist.pdf")
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
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_rminusz_{'density' if density else 'nodensity'}_hist.pdf")
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
plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{REGION}_color.pdf")
#plt.show()
plt.clf()
