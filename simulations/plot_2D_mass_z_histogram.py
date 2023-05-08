import matplotlib
from matplotlib import pyplot as plt
import numpy as np

import directories

from manage_parameter_space import get_details_of_run

plt.rcParams["font.size"] = "16"
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["savefig.pad_inches"] = 0.05
plt.rcParams["savefig.bbox"] = "tight"
#plt.rcParams["figure.figsize"] = (30, 8) # Default: (6.4, 4.8)

pinocchio_particle_count = 2048
Z_DEPTH = 0.5
PINOCCHIO_REGION = "fullsky"

DESI_region = directories.BASS_MzLS

run_ids = [156, 143]

fig, axs = plt.subplots(nrows=len(run_ids), ncols=2, sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.1, hspace=0.15)


for i, (ax_row, run_id) in enumerate(zip(axs, run_ids)):
    # Obtain details of this SHAM run,
    # in particular the desired mass limit
    run_details = get_details_of_run(run_id)
    M_limit = run_details["mass_cut"] # mass limit for assigning blue or red galaxies to halos, [Msun/h]
    num_z_bins = run_details["num_z_bins"]
    num_mass_bins = run_details["num_mass_bins"]

    for ax, color in zip(ax_row, ("blue", "red")):
        # Directory containing output data from SHAM
        infile_dir = directories.path_interpolation(
            particle_count=pinocchio_particle_count,
            z_depth=Z_DEPTH,
            pinocchio_region=PINOCCHIO_REGION,
            DESI_region=DESI_region,
            run_id=run_id,
        )

        infile_hist = f"{color}_lim_interp.npz"

        with np.load(infile_dir + infile_hist) as d:
            lim_abs_mag = d[f"lim_abs_mag_{color}"]
            bin_edges_z = d["bin_edges_z"]
            bin_edges_mass = d["bin_edges_mass"]

        histogram = ax.pcolormesh(bin_edges_z, bin_edges_mass, lim_abs_mag.T, vmin=-24, vmax=-18)

        title_label = f"$m_{{bins}}$ = {num_mass_bins}, $z_{{bins}}$ = {num_z_bins}"
        if i == 0:
            title_label = f"{color.capitalize()} galaxies\n" + title_label
        ax.set_title(title_label)

        ax.set_yscale("log")
        ax.set_ylim([5e10, 1e15])

fig.supxlabel("Halo/Subhalo Redshift")
fig.supylabel(r"Halo/Subhalo mass [M$_{\odot}$/h]")

cbar = fig.colorbar(
    histogram,
    label="Absolute Magnitude",
    ax=axs,
    shrink=0.9,
    fraction=0.046,
    pad=0.04,
)
cbar.ax.invert_yaxis()

fig.savefig(f"/cluster/home/lmachado/msc-thesis/simulations/images/runs{run_ids}_2D_histogram.pdf")

plt.show()
