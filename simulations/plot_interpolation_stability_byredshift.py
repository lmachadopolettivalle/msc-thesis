# For a given parameter space run,
# plot absolute magnitude vs. halo/subhalo mass
# for several different redshifts.
# Make two plots, one for blue and one for red galaxies.

from cycler import cycler
from matplotlib import pyplot as plt
import numpy as np
import re
from scipy.interpolate import griddata
from tqdm import tqdm

import directories

from manage_parameter_space import get_details_of_run

# Colors
blue = "#004488"
yellow = "#ddaa33"
red = "#bb5566"

color_cycler = cycler(color=[
    "#E8601C",
    "#F6C141",
    "#90C987",
    "#5289C7",
    "#AE76A3",
])

plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (12, 7)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"

pinocchio_particle_count = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

DESI_region = directories.BASS_MzLS

# Loop through desired run IDs
DESIRED_RUN_IDS = [156, 143]

# Convert from absolute magnitude to luminosity
def mag_to_luminosity(mag):
    SUN_ABS_MAG = 4.83
    return np.power(10, 0.4 * (SUN_ABS_MAG - mag))

def create_callback(lum_ax):
    def convert_mag_axis_to_luminosity_axis(mag_ax):
        y1, y2 = mag_ax.get_ylim()
        lum_ax.set_ylim(mag_to_luminosity(y1), mag_to_luminosity(y2))
        lum_ax.figure.canvas.draw()

    return convert_mag_axis_to_luminosity_axis

for run_id in tqdm(DESIRED_RUN_IDS):
    # File with output of SHAM run.
    # Used to determine effective mass limit used in SHAM
    sham_output_filename = f"/cluster/home/lmachado/msc-thesis/simulations/sham_int_job_{pinocchio_particle_count}_{run_id}_{DESI_region}_output"

    # Directory containing output data from SHAM
    infile_dir = directories.path_interpolation(
        particle_count=pinocchio_particle_count,
        z_depth=Z_DEPTH,
        pinocchio_region=PINOCCHIO_REGION,
        DESI_region=DESI_region,
        run_id=run_id,
    )
    infile_hist_blue = "blue_lim_interp.npz"
    infile_hist_red = "red_lim_interp.npz"

    # Obtain details of this SHAM run,
    # in particular the desired mass limit
    run_details = get_details_of_run(run_id)
    M_limit = run_details["mass_cut"] # mass limit for assigning blue or red galaxies to halos, [Msun/h]
    num_z_bins = run_details["num_z_bins"]
    num_mass_bins = run_details["num_mass_bins"]

    # Obtain effective mass limit used
    with open(sham_output_filename, 'r') as f:
        text = f.read()

    M_limit_effective = float(
        re.search("M_limit_effective\s+=\s+(.+)\n", text).groups()[0]
    )

    # ----------------------------------------------------–
    # LOAD HALO-SUBHALO HISTOGRAMS
    # -----------------------------------------------------
    with np.load(infile_dir + infile_hist_blue) as data:
        lim_abs_mag_blue = data["lim_abs_mag_blue"]
        blue_bin_edges_z = data["bin_edges_z"]
        blue_bin_edges_mass = data["bin_edges_mass"]

    with np.load(infile_dir + infile_hist_red) as data:
        lim_abs_mag_red = data["lim_abs_mag_red"]
        red_bin_edges_z = data["bin_edges_z"]
        red_bin_edges_mass = data["bin_edges_mass"]

    # Create arrays for interpolation
    blue_mass_edges_stacked = np.reshape(
        np.tile(blue_bin_edges_mass, len(blue_bin_edges_z)),
        (len(blue_bin_edges_mass) * len(blue_bin_edges_z), ),
    )
    red_mass_edges_stacked = np.reshape(
        np.tile(red_bin_edges_mass, len(red_bin_edges_z)),
        (len(red_bin_edges_mass) * len(red_bin_edges_z), ),
    )

    blue_z_edges_stacked = np.reshape(
        np.ndarray.flatten(
            np.transpose(
                np.tile(
                    blue_bin_edges_z,
                    (len(blue_bin_edges_mass), 1),
                )
            )
        ),
        (len(blue_bin_edges_mass)*len(blue_bin_edges_z), ),
    )
    red_z_edges_stacked = np.reshape(
        np.ndarray.flatten(
            np.transpose(
                np.tile(
                    red_bin_edges_z,
                    (len(red_bin_edges_mass), 1),
                )
            )
        ),
        (len(red_bin_edges_mass)*len(red_bin_edges_z), ),
    )

    # Choose input values for mass and redshift
    input_masses = np.linspace(11, 15, num=100)
    input_redshifts = np.linspace(0.05, 0.5, 5)
    
    fig, (mag_ax_blue, mag_ax_red) = plt.subplots(nrows=1, ncols=2, sharey=True)
    fig.subplots_adjust(wspace=0.05)
    mag_ax_blue.set_prop_cycle(color_cycler)
    mag_ax_red.set_prop_cycle(color_cycler)
    lum_ax_blue = mag_ax_blue.twinx()
    lum_ax_red = mag_ax_red.twinx()

    # automatically update ylim of ax2 when ylim of ax1 changes.
    mag_ax_blue.callbacks.connect("ylim_changed", create_callback(lum_ax_blue))
    mag_ax_red.callbacks.connect("ylim_changed", create_callback(lum_ax_red))

    lum_ax_blue.set_yscale("log")
    lum_ax_red.set_yscale("log")

    # Interpolate for each redshift
    for z in input_redshifts:
        blue_absmags = griddata(
            (blue_z_edges_stacked, np.log10(blue_mass_edges_stacked)),
            np.ndarray.flatten(lim_abs_mag_blue),
            (z, input_masses),
        )
        red_absmags = griddata(
            (red_z_edges_stacked, np.log10(red_mass_edges_stacked)),
            np.ndarray.flatten(lim_abs_mag_red),
            (z, input_masses),
        )

        mag_ax_blue.plot(
            10**input_masses,
            blue_absmags,
            label=f"z = {z:.2f}",
            linewidth=2,
        )
        mag_ax_red.plot(
            10**input_masses,
            red_absmags,
            label=f"z = {z:.2f}",
            linewidth=2,
        )

    for color, ax in zip(["blue", "red"], [mag_ax_blue, mag_ax_red]):
        ax.invert_yaxis()

        ax.set_xscale("log")

        ax.set_ylim([-18, -24])

        ax.set_xlim([5e10, 1e15])


        ax.text(
            1e12,
            -23.5,
            f"{color.capitalize()} Galaxies\n$m_{{bins}}$ = {num_mass_bins}, $z_{{bins}}$ = {num_z_bins}",
            color=(blue if color == "blue" else red),
            bbox=dict(facecolor="none", edgecolor="black", pad=6),
            ha="center",
            va="center",
        )

        ax.grid()

    fig.supxlabel(r"Halo/Subhalo Mass ($M_{\odot}/h$)")
    mag_ax_blue.set_ylabel("Galaxy Absolute Magnitude")

    lum_ax_red.set_ylabel("Galaxy Luminosity ($L_{\odot}$)")

    mag_ax_blue.legend(loc="lower right")
    lum_ax_blue.yaxis.set_ticks_position("none") 
    lum_ax_blue.tick_params(
        axis="y",          # changes apply to the x-axis
        which="both",      # both major and minor ticks are affected
        left=False,      # ticks along the left edge are off
        right=False,         # ticks along the right edge are off
        labelright=False,
    )

    fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/interpolation_byredshift_{pinocchio_particle_count}_{run_id}.pdf")

    plt.show()

    # Close figures to reduce memory usage,
    # since pyplot keeps figures in memory until
    # the end of the program by default
    plt.close(fig)
