from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate

import directories
from manage_parameter_space import get_details_of_run

nbins = 7
TMP_BINS = np.logspace(np.log10(0.06), np.log10(3), nbins)
BINS = (TMP_BINS[:-1] + TMP_BINS[1:]) / 2
BINS = BINS[:-1]

# Number of particles (cube root) used in run
# This determines the path where the data is stored
PARTICLE_COUNT_PINOCCHIO = 2048
Z_DEPTH = 0.5

RUN_ID = 148

DESI_REGION = directories.FULLSKY
PINOCCHIO_REGION = "fullsky"

if DESI_REGION == directories.FULLSKY:
    EQUIVALENT_DESI_REGION = (directories.BASS_MzLS, directories.DECaLS_NGC, directories.DECaLS_SGC)
else:
    EQUIVALENT_DESI_REGION = DESI_REGION

PATH_DESI_LS_2PCF = "/cluster/scratch/lmachado/DataProducts/2PCF/"

rmag_bins = [
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 19.5],
]

# Angular scales of trust for each r band range
SCALES_OF_TRUST = {
    (15, 16): (0.26, 3),
    (16, 17): (0.20, 3),
    (17, 18): (0.1, 3),
    (18, 19): (0.09, 3),
    (19, 19.5): (0.08, 3),
}

measurements = None

for seed in (
    None,
    "111111",
    "222222",
    "333333",
    "444444",
    "555555",
    "666666",
):
    PATH_SHAM_2PCF = directories.path_2PCF(
        particle_count=PARTICLE_COUNT_PINOCCHIO,
        z_depth=Z_DEPTH,
        pinocchio_region=PINOCCHIO_REGION,
        DESI_region=DESI_REGION,
        run_id=RUN_ID,
        seed=seed,
    )

    rmag_bin_measurements = None

    for rmag_low, rmag_high in rmag_bins:
        sham_bins_filename = f"simulated_total_2PCF_rprimed_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_rprimed_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
            sham_bins = np.load(f)
        with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
            sham_wtheta = np.load(f)

        # Only use bins for which we trust the 2PCF
        trusted_indices = np.where(
            (sham_bins >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
            (sham_bins <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
        )[0]
        BINS_trusted = BINS[
            (BINS >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
            (BINS <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
        ]

        sham_wtheta = interpolate.interp1d(
            sham_bins[trusted_indices],
            sham_wtheta[trusted_indices],
        )(
            BINS_trusted
        )
        sham_bins = BINS_trusted
        plt.plot(sham_bins, sham_wtheta, label=f"{rmag_low:.1f} - {rmag_high:.1f}")

        if rmag_bin_measurements is None:
            rmag_bin_measurements = sham_wtheta
        else:
            rmag_bin_measurements = np.hstack([rmag_bin_measurements, sham_wtheta])

    if measurements is None:
        measurements = rmag_bin_measurements
    else:
        measurements = np.vstack([measurements, rmag_bin_measurements])

# Compute covariance matrix
mean = np.mean(measurements, axis=0)
cov = np.cov((measurements - mean).T)
#cov_inverse = np.linalg.inv(cov) # Unstable inversion, should obtain via regression instead

#plt.matshow(cov)
#plt.matshow(cov_inverse)
#plt.show()

# Load DESI LS 2PCF for chi square computation
data = None
for rmag_low, rmag_high in rmag_bins:
    desi_bins_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_bins_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_primed.npy"
    desi_wtheta_filename = f"{EQUIVALENT_DESI_REGION}_2PCF_wtheta_Bright_rmag_range{rmag_low:.1f}-{rmag_high:.1f}_primed.npy"

    with open(f"{PATH_DESI_LS_2PCF}/{desi_bins_filename}", "rb") as f:
        desi_bins = np.load(f)
    with open(f"{PATH_DESI_LS_2PCF}/{desi_wtheta_filename}", "rb") as f:
        desi_wtheta = np.load(f)

    trusted_indices = np.where(
        (desi_bins >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
        (desi_bins <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
    )[0]
    BINS_trusted = BINS[
        (BINS >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
        (BINS <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
    ]

    desi_wtheta = interpolate.interp1d(
        desi_bins[trusted_indices],
        desi_wtheta[trusted_indices],
        fill_value="extrapolate",
    )(
        BINS_trusted
    )
    desi_bins = BINS_trusted
    plt.plot(desi_bins, 1.06*desi_wtheta, ls="--", lw=5, label="DATA")


    if data is None:
        data = desi_wtheta
    else:
        data = np.hstack((data, desi_wtheta))

data *= 1.06

plt.legend()
plt.show()
# Difference between simulation and data,
# used for chi square
for m in measurements:
    difference = data - m

    # Use regression to estimate inverse of covariance matrix and compute chi square
    lambda_param = 1e-18
    chi_square_approx = np.dot(
        difference.T,
        np.linalg.solve(
            np.dot(cov.T, cov) + lambda_param * np.identity(len(difference)),
            np.dot(cov.T, difference)
        )
    )
    print(chi_square_approx)

    # Compare against chi square obtained from directly inverting covariance matrix
    #chi_square_direct = np.matmul(difference, np.matmul(cov_inverse, difference.T))
    #print(chi_square_direct)

print("Computing chi square for all SHAM parameter values")

chi_square_dict = {} # Key = run_id, Value = chi square
run_details_dict = {} # Key = run_id, Value = run details

RUN_IDS = [140, 141, 142, 143] + list(range(146, 155 + 1)) + list(range(157, 172 + 1))

for run_id in RUN_IDS:
    PATH_SHAM_2PCF = directories.path_2PCF(
        particle_count=PARTICLE_COUNT_PINOCCHIO,
        z_depth=Z_DEPTH,
        pinocchio_region=PINOCCHIO_REGION,
        DESI_region=DESI_REGION,
        run_id=run_id,
        seed=None,
    )

    rmag_bin_measurements = None

    for rmag_low, rmag_high in rmag_bins:
        sham_bins_filename = f"simulated_total_2PCF_rprimed_{rmag_low:.1f}_{rmag_high:.1f}_bins.npy"
        sham_wtheta_filename = f"simulated_total_2PCF_rprimed_{rmag_low:.1f}_{rmag_high:.1f}_wtheta.npy"

        with open(f"{PATH_SHAM_2PCF}/{sham_bins_filename}", "rb") as f:
            sham_bins = np.load(f)
        with open(f"{PATH_SHAM_2PCF}/{sham_wtheta_filename}", "rb") as f:
            sham_wtheta = np.load(f)

        # Only use bins for which we trust the 2PCF
        trusted_indices = np.where(
            (sham_bins >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
            (sham_bins <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
        )[0]
        BINS_trusted = BINS[
            (BINS >= SCALES_OF_TRUST[(rmag_low, rmag_high)][0]) &
            (BINS <= SCALES_OF_TRUST[(rmag_low, rmag_high)][1])
        ]

        sham_wtheta = interpolate.interp1d(
            sham_bins[trusted_indices],
            sham_wtheta[trusted_indices],
            fill_value="extrapolate",
        )(
            BINS_trusted
        )
        sham_bins = BINS_trusted

        if rmag_bin_measurements is None:
            rmag_bin_measurements = sham_wtheta
        else:
            rmag_bin_measurements = np.hstack([rmag_bin_measurements, sham_wtheta])

    difference = data - rmag_bin_measurements
    chi_square_approx = np.dot(
        difference.T,
        np.linalg.solve(
            np.dot(cov.T, cov) + lambda_param * np.identity(len(difference)),
            np.dot(cov.T, difference)
        )
    )

    run_details_dict[run_id] = get_details_of_run(run_id)
    print(f"{run_id}: {run_details_dict[run_id]}. CHI: {chi_square_approx}")

    chi_square_dict[run_id] = chi_square_approx

# Plot chi square values to find minimal configuration
masses = [run_details_dict[run_id]["mass_cut"] for run_id in RUN_IDS]
t_quench = [run_details_dict[run_id]["quenching_time"] for run_id in RUN_IDS]
chi_squares = [chi_square_dict[run_id] for run_id in RUN_IDS]
sc = plt.scatter(
    masses,
    t_quench,
    c=chi_squares,
)
plt.xlabel(r"$M_{limit}$ ($M_{\odot}$)", fontsize=15)
plt.ylabel(r"$t_{quench}$ (Gyr)", fontsize=15)
plt.xscale("log")
plt.colorbar(sc, label=r"$\chi^2$")
plt.show()

# Make interpolated grid of chi square values
from matplotlib import tri
xi = np.logspace(np.log10(1e12), np.log10(2e13), 50)
yi = np.linspace(0, 3, 50)

triang = tri.Triangulation(masses, t_quench)
interpolator = tri.LinearTriInterpolator(triang, chi_squares)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)

plt.contour(xi, yi, zi, levels=6, linewidths=0.5, colors='k')
contour = plt.contourf(xi, yi, zi, levels=6, cmap="RdBu_r")

plt.xscale("log")

plt.colorbar(contour, label=r"$\chi^2$")
plt.xlabel(r"$M_{limit}$ ($M_{\odot}$)", fontsize=15)
plt.ylabel(r"$t_{quench}$ (Gyr)", fontsize=15)

plt.show()
