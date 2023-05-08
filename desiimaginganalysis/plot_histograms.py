# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np

from constants import BANDS as CONSTANT_BANDS, BGS_BRIGHT, BGS_FAINT, BASS_MzLS, DECaLS_NGC, DECaLS_SGC, map_region_to_north_south
from load_processed_target_data import load_processed_target_data

# Make copy of CONSTANT_BANDS to allow for local modifications
BANDS = CONSTANT_BANDS.copy()

# Whether to plot faint targets
PLOT_FAINT_TARGETS = False

# Add r_primed as a band for histogram plots
# Note that we use the unprimed r for color computations,
# since there is no primed g or z
BANDS.append("r_primed")
BANDS.sort()

# If True, use magnitudes with extinction correction
# If False, show magnitudes without extinction correction
extinction_correction = True

# Parameters for plotting
bin_count = 90

LINEWIDTH = 1.5

ALPHA = 0.5

# Color choices
blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"
orange = "#EE7733"
teal = "#009988"

plt.rcParams["font.size"] = "16"
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"

bright_plot_color = blue
faint_plot_color = orange

HISTTYPE = {
    "north": "step",
    "south": "stepfilled",
}

# Number of pixels in each region,
# to help normalize histograms
NUMBER_PIXELS_AFTER_MASKING = {
    "north": 24684,
    "south": 68236,
}

# Load targets
print("Loading target data")
targets = {
    "north": load_processed_target_data(regions={BASS_MzLS}, extinction_correction=extinction_correction, apply_mask=True),
    "south": load_processed_target_data(regions={DECaLS_NGC, DECaLS_SGC}, extinction_correction=extinction_correction, apply_mask=True),
}

bright_targets = {
    region: np.where(values["BGS_TARGET"] == BGS_BRIGHT)[0]
    for region, values in targets.items()
}
faint_targets = {
    region: np.where(values["BGS_TARGET"] == BGS_FAINT)[0]
    for region, values in targets.items()
}

# Plot histograms of magnitudes in 3 bands
print("Number of north bright objects:", len(targets["north"]["MAG_R"][bright_targets["north"]]))
print("Number of north faint objects:", len(targets["north"]["MAG_R"][faint_targets["north"]]))
print("Number of south bright objects:", len(targets["south"]["MAG_R"][bright_targets["south"]]))
print("Number of south faint objects:", len(targets["south"]["MAG_R"][faint_targets["south"]]))

# Determine binning for each color band
bins = {
    "g": np.linspace(15, 22, bin_count),
    "r_primed": np.linspace(14, 20.5, bin_count),
    "r": np.linspace(14, 20.5, bin_count),
    "z": np.linspace(15, 20.5, bin_count),
}

# Begin plotting
print("Plotting magnitude histograms")

# Overplot r-primed and r-unprimed for comparison
bright_label = "BGS Bright"
faint_label = "BGS Faint"

values = targets["north"]

bright_r_values = values["MAG_R"][bright_targets["north"]]
bright_r_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(bright_r_values)

bright_rprimed_values = values["MAG_R_PRIMED"][bright_targets["north"]]
bright_rprimed_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(bright_rprimed_values)

faint_r_values = values["MAG_R"][faint_targets["north"]]
faint_r_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(faint_r_values)

faint_rprimed_values = values["MAG_R_PRIMED"][faint_targets["north"]]
faint_rprimed_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(faint_rprimed_values)

plt.hist(
    bright_r_values,
    weights=bright_r_weights,
    alpha=1,
    bins=bins["r"], linewidth=LINEWIDTH, density=False, histtype="step", label=f"{bright_label}, unprimed", color=bright_plot_color
)
if PLOT_FAINT_TARGETS:
    plt.hist(
        faint_r_values,
        weights=faint_r_weights,
        alpha=1,
        bins=bins["r"], linewidth=LINEWIDTH, density=False, histtype="step", label=f"{faint_label}, unprimed", color=faint_plot_color
    )

plt.hist(
    bright_rprimed_values,
    weights=bright_rprimed_weights,
    alpha=ALPHA,
    bins=bins["r_primed"], linewidth=LINEWIDTH, density=False, histtype="stepfilled", label=f"{bright_label}, primed", color=bright_plot_color
)
if PLOT_FAINT_TARGETS:
    plt.hist(
        faint_rprimed_values,
        weights=faint_rprimed_weights,
        alpha=ALPHA,
        bins=bins["r_primed"], linewidth=LINEWIDTH, density=False, histtype="stepfilled", label=f"{faint_label}, primed", color=faint_plot_color
    )

plt.xlim([14, 22])

#plt.title("Comparison r-primed and unprimed, BASS")
plt.xlabel("r magnitude, extinction corrected")

plt.ylabel("Count, normalized by survey area")
plt.legend(loc="upper left")
plt.grid()
plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/r_primed_r_comparison.pdf")

plt.show()
#plt.clf()

# Plot histograms
fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(12.5, 5))
fig.subplots_adjust(wspace=0.15)
for ax, band in zip(axs, ["g", "r_primed", "z"]):
    print(band)

    if extinction_correction:
        bright_label = "BGS Bright"
        faint_label = "BGS Faint"
    else:
        bright_label = "BGS Bright (E.C.)"
        faint_label = "BGS Faint (E.C.)"

    for region, values in targets.items():
        bright_values = values[f"MAG_{band.upper()}"][bright_targets[region]]
        bright_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(bright_values)

        faint_values = values[f"MAG_{band.upper()}"][faint_targets[region]]
        faint_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(faint_values)

        ax.hist(
            bright_values,
            #weights=bright_weights,
            alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
            bins=bins[band], linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"{bright_label}, {region}", color=bright_plot_color
        )
        if PLOT_FAINT_TARGETS:
            ax.hist(
                faint_values,
                #weights=faint_weights,
                alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
                bins=bins[band], linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"{faint_label}, {region}", color=faint_plot_color
            )

    ax.set_xlim([14, 22])
    if band == "g":
        ax.legend(loc="upper left")
    ax.grid()

    ax.set_xlabel(f"{band[0]} mag")

fig.supylabel("PDF")

if extinction_correction:
    fig.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/mag_hist.pdf")
else:
    fig.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/mag_hist_no_extinction.pdf")

plt.show()




# Plot histograms of colors
bright_gminusr = {
    region: values["MAG_G"][bright_targets[region]] - values["MAG_R"][bright_targets[region]]
    for region, values in targets.items()
}
faint_gminusr = {
    region: values["MAG_G"][faint_targets[region]] - values["MAG_R"][faint_targets[region]]
    for region, values in targets.items()
}

bright_rminusz = {
    region: values["MAG_R"][bright_targets[region]] - values["MAG_Z"][bright_targets[region]]
    for region, values in targets.items()
}
faint_rminusz = {
    region: values["MAG_R"][faint_targets[region]] - values["MAG_Z"][faint_targets[region]]
    for region, values in targets.items()
}

### Color histograms
bins_color_histograms = np.linspace(-1, 2, bin_count)
# G - R Color Histogram
fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8, 5))
fig.subplots_adjust(wspace=0.15)
print("Plotting G - R color histogram")
for region in bright_gminusr.keys():
    bright_values = bright_gminusr[region]
    bright_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(bright_values)
    faint_values = faint_gminusr[region]
    faint_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(faint_values)

    axs[0].hist(
        bright_values,
        #weights=bright_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"BGS Bright, {region}", color=bright_plot_color
    )
    if PLOT_FAINT_TARGETS:
        axs[0].hist(
            faint_values,
            #weights=faint_weights,
            alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
            bins=bins_color_histograms, linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"BGS Faint, {region}", color=faint_plot_color
        )

axs[0].set_xlim([0, 2])
axs[0].set_xlabel("g - r")
axs[0].grid()

# R - Z Color Histogram
print("Plotting R - Z color histogram")
for region in bright_rminusz.keys():
    bright_values = bright_rminusz[region]
    bright_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(bright_values)
    faint_values = faint_rminusz[region]
    faint_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(faint_values)

    axs[1].hist(
        bright_values,
        #weights=bright_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"BGS Bright, {region}", color=bright_plot_color
    )
    if PLOT_FAINT_TARGETS:
        axs[1].hist(
            faint_values,
            #weights=faint_weights,
            alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
            bins=bins_color_histograms, linewidth=LINEWIDTH, density=True, histtype=HISTTYPE[region], label=f"BGS Faint, {region}", color=faint_plot_color
        )

axs[1].set_xlim([0, 2])
axs[1].set_xlabel("r - z")
axs[1].grid()

fig.supylabel("PDF")

fig.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/colors_hist.pdf")

plt.show()

exit()

#####
print("Plotting color scatter plot")
# G-R vs. R-Z color plot
# Since this is a scatter plot,
# cannot handle plotting all targets.
# Instead, sample a subset of the targets to be plotted.
sample_number = 10000
for region in bright_rminusz.keys():
    bright_sample_idx = np.random.choice(np.arange(len(bright_rminusz[region])), sample_number, replace=False)
    faint_sample_idx = np.random.choice(np.arange(len(faint_rminusz[region])), sample_number, replace=False)

    plt.scatter(bright_rminusz[region][bright_sample_idx], bright_gminusr[region][bright_sample_idx], label=f"BGS Bright", s=4, alpha=0.5, c=bright_plot_color)
    if PLOT_FAINT_TARGETS:
        plt.scatter(faint_rminusz[region][faint_sample_idx], faint_gminusr[region][faint_sample_idx], label=f"BGS Faint", s=4, alpha=0.5, c=faint_plot_color)

    plt.xlabel("r - z color")
    plt.ylabel("g - r color")
    plt.xlim([-1.5, 3.5])
    plt.ylim([-1.5, 3.5])
    #plt.title(f"Color plot for {region} region")
    plt.legend(loc="upper right")
    plt.grid()
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/color_{region}_unprimed.pdf")
    #plt.show()
    plt.clf()
