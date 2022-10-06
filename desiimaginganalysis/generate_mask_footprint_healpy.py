from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
from shapely.geometry import LineString, Point, Polygon

from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned

from generate_mask_footprint_shapely import load_boundary_multipolygon

# Load shapely footprint polygons
shapely_footprint = load_boundary_multipolygon()

# Create healpy pixels with a given nside
nside = 4
npix = hp.nside2npix(nside)
print("npix: ", npix)
pixel_indices = range(npix)

print(np.degrees(hp.nside2resol(nside)))

pixels_in_footprint = []

# Loop through the pixels, and get the coordinates of their boundary points
for pixel_id in pixel_indices:
    # Obtain boundary coords in (x, y, z)
    pixel_boundary_coords = np.transpose(
        hp.boundaries(
            nside,
            pixel_id,
            step=1, # Only retrieve the 4 corners
            nest=False # Decide between RING and NEST
        )
    )
    # Convert to RA and Dec (RA: 0-360, DecL -90 - +90)
    ras, decs = hp.pixelfunc.vec2ang(pixel_boundary_coords, lonlat=True)
    print(ras, decs)

    # TODO handle case when the pixel crosses the 0/360 RA line
    # NOTE I believe that, by checking the corners only, this issue is resolved

    points = [Point(r, d) for r, d in zip(ras, decs)]
    contains = True
    for p in points:
        if not shapely_footprint.contains(p):
            contains = False
            break
    if contains:
        pixels_in_footprint.append(pixel_id)


# Save pixelids to file

with open(f"pixels_in_footprint_for_nside_{nside}.npy", "wb") as f:
    np.save(f, np.array(pixels_in_footprint))


ax = init_sky(galactic_plane_color=None, ecliptic_plane_color=None)

m = np.zeros(npix)
m[pixels_in_footprint] = 10 # Some non-zero value
plot_healpix_map(m, nest=False, cmap='viridis', colorbar=True, label=None, ax=ax)

plt.show()
