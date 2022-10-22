import healpy as hp
import numpy as np

RNG = np.random.default_rng(42)

def draw_randoms(number_of_points=100, nside=128):

    U, V = np.array(list(zip(*RNG.random((number_of_points, 2)))))

    RA = np.degrees(2 * np.pi * U)
    DEC = np.degrees(np.pi / 2 - np.arccos(2 * V - 1))


    with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_RA.npy", "wb") as f:
        np.save(f, RA)
    with open("/cluster/scratch/lmachado/DataProducts/randoms/randoms_DEC.npy", "wb") as f:
        np.save(f, DEC)

    pixel_ids = hp.ang2pix(nside, RA, DEC, nest=True, lonlat=True)
    with open(f"/cluster/scratch/lmachado/DataProducts/randoms/randoms_pixels_NSIDE_{nside}.npy", "wb") as f:
        np.save(f, pixel_ids)


if __name__ == "__main__":
    draw_randoms(2 * 30000000, nside=256)
