from matplotlib import pyplot as plt
import numpy as np

from constants import BANDS, BGS_BRIGHT, BGS_FAINT, BASS_MzLS, DECaLS_NGC, DECaLS_SGC
from load_processed_target_data import load_processed_target_data

# Parameters for plotting
ALPHA = 0.6

plt.rcParams['font.size'] = '12'

UNIQUE_MORPHTYPES = ["PSF", "DEV", "EXP", "REX", "SER"]

# Load targets
print("Loading target data")
targets = {
    "north": load_processed_target_data(regions={BASS_MzLS}, extinction_correction=True, apply_mask=True),
    "south": load_processed_target_data(regions={DECaLS_NGC, DECaLS_SGC}, extinction_correction=True, apply_mask=True),
}

for region in targets.keys():
    targets[region]["GMINUSR"] = targets[region]["MAG_G"] - targets[region]["MAG_R"]
    targets[region]["RMINUSZ"] = targets[region]["MAG_R"] - targets[region]["MAG_Z"]


bright_targets = {
    region: np.where(values["BGS_TARGET"] == BGS_BRIGHT)[0]
    for region, values in targets.items()
}
faint_targets = {
    region: np.where(values["BGS_TARGET"] == BGS_FAINT)[0]
    for region, values in targets.items()
}
morphtype_targets = {
    region: {
        morphtype: np.where(values["MORPHTYPE"] == morphtype)[0]
        for morphtype in UNIQUE_MORPHTYPES
    }
    for region, values in targets.items()
}

#####
print("Plotting color scatter plot for each target morphtype")

for region in targets.keys():
    for brightness, brightness_ids in (("bright", bright_targets[region]), ("faint", faint_targets[region])):
        for morphtype in UNIQUE_MORPHTYPES:
            mask = np.intersect1d(
                brightness_ids,
                morphtype_targets[region][morphtype],
            )

            sample_number = 2000
            sample_target_ids = np.random.choice(
                mask,
                sample_number,
                replace=False,
            )

            plt.scatter(
                targets[region]["RMINUSZ"][sample_target_ids],
                targets[region]["GMINUSR"][sample_target_ids],
                label=morphtype,
                s=4,
                alpha=ALPHA,

            )

        plt.xlabel("r - z color")
        plt.ylabel("g - r color")
        plt.xlim([-1.5, 3.5])
        plt.ylim([-1.5, 3.5])
        plt.title(f"Color plot for {region} region, BGS {brightness}")
        plt.legend(loc="upper right")
        plt.grid()
        plt.savefig(f"/cluster/home/lmachado/msc-thesis/desiimaginganalysis/images/morphtypes_color_{brightness}_{region}.pdf")
        plt.show()
