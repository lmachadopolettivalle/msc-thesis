from matplotlib import pyplot as plt
import numpy as np

# Parameters for plotting
plt.rcParams["font.size"] = "12"
LINEWIDTH = 2

rmag_bins = [
    [14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

for color_name in ("blue", "red"):
    for rmag_low, rmag_high in rmag_bins:
        with open(f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy", "rb") as f:
            bins = np.load(f)
        with open(f"simulated_{color_name}_2PCF_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy", "rb") as f:
            wtheta = np.load(f)

        plt.plot(
            bins,
            wtheta,
            linewidth=LINEWIDTH,
            label=f"{rmag_low:.1f} < r < {rmag_high:.1f}",
        )

    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\theta$ [deg]")
    plt.ylabel(r"$w(\theta)$")
    plt.title(f"Simulated galaxies, {color_name}")
    plt.xlim([2e-3, 20])
    plt.ylim([1e-4, 100])
    plt.grid()

    #plt.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/2PCF_BASS_simulated_{color_name}.pdf")

    plt.show()
