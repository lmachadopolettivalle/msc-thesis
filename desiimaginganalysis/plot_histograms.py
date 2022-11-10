# Read file with data from selected targets
# Plot histogram of magnitudes and colors
# Filter based on BGS bright vs. BGS faint

from matplotlib import pyplot as plt
import numpy as np

from constants import BANDS as CONSTANT_BANDS, BGS_BRIGHT, BGS_FAINT, BASS_MzLS, DECaLS_NGC, DECaLS_SGC, map_region_to_north_south
from load_processed_target_data import load_processed_target_data

# Make copy of CONSTANT_BANDS to allow for local modifications
BANDS = CONSTANT_BANDS.copy()

# Choose whether to use the primed or unprimed version for the r-mag values
MAG_R_PRIMED = False

if MAG_R_PRIMED:
    BANDS.remove("r")
    BANDS.append("r_primed")
    BANDS.sort()

if MAG_R_PRIMED:
    MAG_R = "MAG_R_PRIMED"
else:
    MAG_R = "MAG_R"

# If True, use magnitudes with extinction correction
# If False, show magnitudes without extinction correction
extinction_correction = True

# Parameters for plotting
bin_count = 90

LINEWIDTH = 1.25

ALPHA = 0.5

# Color choices
blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"

plt.rcParams['font.size'] = '12'
bright_plot_color = blue
faint_plot_color = red

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
print("Number of north bright objects:", len(targets["north"][MAG_R][bright_targets["north"]]))
print("Number of north faint objects:", len(targets["north"][MAG_R][faint_targets["north"]]))
print("Number of south bright objects:", len(targets["south"][MAG_R][bright_targets["south"]]))
print("Number of south faint objects:", len(targets["south"][MAG_R][faint_targets["south"]]))

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
bright_label = "BGS Bright (E.C.)"
faint_label = "BGS Faint (E.C.)"

values = targets["north"]

bright_r_values = values["MAG_R"][bright_targets["north"]]
bright_r_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(bright_values)

bright_rprimed_values = values["MAG_R_PRIMED"][bright_targets["north"]]
bright_rprimed_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(bright_values)

faint_r_values = values["MAG_R"][faint_targets["north"]]
faint_r_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(faint_values)

faint_rprimed_values = values["MAG_R_PRIMED"][faint_targets["north"]]
faint_rprimed_weights = [1 / NUMBER_PIXELS_AFTER_MASKING["north"]] * len(faint_values)

plt.hist(
    bright_r_values,
    weights=bright_r_weights,
    alpha=1,
    bins=bins["r"], linewidth=LINEWIDTH, density=False, histtype="step", label=f"{bright_label}, unprimed", color=bright_plot_color
)
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
plt.hist(
    faint_rprimed_values,
    weights=faint_rprimed_weights,
    alpha=ALPHA,
    bins=bins["r_primed"], linewidth=LINEWIDTH, density=False, histtype="stepfilled", label=f"{faint_label}, primed", color=faint_plot_color
)

plt.xlim([14, 22])

plt.title("Comparison r-primed and unprimed, BASS")
plt.xlabel("r magnitude, extinction corrected")

plt.ylabel("Count, normalized by survey area")
plt.legend(loc="upper left")
plt.grid()
plt.savefig("/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/r_primed_r_comparison.pdf")

#plt.show()
plt.clf()

exit() # TODO remove me


# Plot histograms
for band in BANDS:
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

        plt.hist(
            bright_values,
            weights=bright_weights,
            alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
            bins=bins[band], linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"{bright_label}, {region}", color=bright_plot_color
        )
        plt.hist(
            faint_values,
            weights=faint_weights,
            alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
            bins=bins[band], linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"{faint_label}, {region}", color=faint_plot_color
        )

    plt.xlim([14, 22])

    if extinction_correction:
        plt.title(f"{band} mag distribution")
        plt.xlabel(f"{band} magnitude, extinction-corrected")
    else:
        plt.title(f"{band} mag, with and without extinction correction (E.C.)")
        plt.xlabel(f"{band} magnitude")

    plt.ylabel("Count, normalized by survey area")
    plt.legend(loc="upper left")
    plt.grid()
    if extinction_correction:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{band}mag_hist.pdf")
    else:
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/{band}mag_hist_no_extinction.pdf")

    #plt.show()
    plt.clf()

# Plot histograms of colors
bright_gminusr = {
    region: values["MAG_G"][bright_targets[region]] - values[MAG_R][bright_targets[region]]
    for region, values in targets.items()
}
faint_gminusr = {
    region: values["MAG_G"][faint_targets[region]] - values[MAG_R][faint_targets[region]]
    for region, values in targets.items()
}

bright_rminusz = {
    region: values[MAG_R][bright_targets[region]] - values["MAG_Z"][bright_targets[region]]
    for region, values in targets.items()
}
faint_rminusz = {
    region: values[MAG_R][faint_targets[region]] - values["MAG_Z"][faint_targets[region]]
    for region, values in targets.items()
}

### Color histograms
bins_color_histograms = np.linspace(-1, 2, bin_count)
# G - R Color Histogram
print("Plotting G - R color histogram")
for region in bright_gminusr.keys():
    bright_values = bright_gminusr[region]
    bright_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(bright_values)
    faint_values = faint_gminusr[region]
    faint_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(faint_values)

    plt.hist(
        bright_values,
        weights=bright_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"BGS Bright, {region}", color=bright_plot_color
    )
    plt.hist(
        faint_values,
        weights=faint_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"BGS Faint, {region}", color=faint_plot_color
    )

plt.xlim([0, 2])
plt.xlabel("g - r color")
plt.ylabel("Count, normalized by survey area")
plt.legend(loc="upper right")
plt.grid()
plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/gminusr_{'primed' if MAG_R_PRIMED else 'unprimed'}_hist.pdf")
#plt.show()
plt.clf()

# R - Z Color Histogram
print("Plotting R - Z color histogram")
for region in bright_rminusz.keys():
    bright_values = bright_rminusz[region]
    bright_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(bright_values)
    faint_values = faint_rminusz[region]
    faint_weights = [1 / NUMBER_PIXELS_AFTER_MASKING[region]] * len(faint_values)

    plt.hist(
        bright_values,
        weights=bright_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"BGS Bright, {region}", color=bright_plot_color
    )
    plt.hist(
        faint_values,
        weights=faint_weights,
        alpha=(ALPHA if "filled" in HISTTYPE[region] else 1),
        bins=bins_color_histograms, linewidth=LINEWIDTH, density=False, histtype=HISTTYPE[region], label=f"BGS Faint, {region}", color=faint_plot_color
    )

plt.xlim([0, 2])
plt.xlabel("r - z color")
plt.ylabel("Count, normalized by survey area")
plt.legend(loc="upper right")
plt.grid()
plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/rminusz_{'primed' if MAG_R_PRIMED else 'unprimed'}_hist.pdf")
#plt.show()
plt.clf()

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
    plt.scatter(faint_rminusz[region][faint_sample_idx], faint_gminusr[region][faint_sample_idx], label=f"BGS Faint", s=4, alpha=0.5, c=faint_plot_color)

    plt.xlabel("r - z color")
    plt.ylabel("g - r color")
    plt.xlim([-1.5, 3.5])
    plt.ylim([-1.5, 3.5])
    plt.title(f"Color plot for {region} region")
    plt.legend(loc="upper right")
    plt.grid()
    plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/color_{region}_{'primed' if MAG_R_PRIMED else 'unprimed'}.pdf")
    #plt.show()
    plt.clf()
