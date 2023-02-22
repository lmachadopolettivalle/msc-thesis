import numpy as np

import directories

# NOTE: very careful when loading both coordinates (x, y, z)
# and bands (g, r, z), since there is an overlap of names (z).
# Make sure to name them differently, e.g. mag_z and z_coord
BANDS = ["mag_g", "mag_r", "mag_z"]

def load_sham_galaxies(
    particle_count=2048,
    z_depth=0.5,
    pinocchio_region="fullsky",
    DESI_region=directories.BASS_MzLS,
    run_id=146,
):
    # Load x, y, z positions
    # Path to output data from SHAM
    SHAM_OUTPUT_PATH = directories.path_interpolation(
        particle_count=particle_count,
        z_depth=z_depth,
        pinocchio_region=pinocchio_region,
        DESI_region=DESI_region,
        run_id=run_id,
    )

    galaxies = {}
    for coord in ("x_coord", "y_coord", "z_coord"):
        filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_{coord}.npy"

        galaxies[coord] = np.load(filename)

    # Load magnitudes
    for band in BANDS:
        filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_app_{band}.npy"

        galaxies[band] = np.load(filename)

    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_z.npy"
    galaxies["redshift"] = np.load(filename)

    # Load whether galaxies are blue or red
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_blue_red.npy"

    galaxies["blue_red"] = np.load(filename)

    # Load absolute magnitudes
    filename = f"{SHAM_OUTPUT_PATH}/ucat_sorted_app_mag_interp_abs_mag.npy"

    galaxies["abs_mag"] = np.load(filename)

    # Convert 3D positions into RA, Dec
    radii = np.sqrt(galaxies["x_coord"]**2 + galaxies["y_coord"]**2 + galaxies["z_coord"]**2)

    theta = np.arccos(galaxies["z_coord"] / radii)
    phi = np.arctan2(galaxies["y_coord"], galaxies["x_coord"])

    # Note that phi is in range [-pi, pi], but for healpy, must be in range [0, 360 degrees]
    phi[phi < 0] += 2 * np.pi

    galaxies["RA"] = np.degrees(phi)
    galaxies["DEC"] = np.degrees(np.pi/2 - theta)
    
    return galaxies
