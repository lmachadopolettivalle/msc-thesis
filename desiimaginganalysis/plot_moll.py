# Guide: https://github.com/desihub/desiutil/blob/main/doc/nb/SkyMapExamples.ipynb
# Docs: https://desiutil.readthedocs.io/en/latest/api.html?highlight=init_sky#desiutil.plots.init_sky
# See usage examples in QA file: https://github.com/desihub/desitarget/blob/34e5f0322a96803896392d32131a685d007b12cc/py/desitarget/QA.py
# For help with mask:
# See also _in_desi_footprint

from desiutil.plots import prepare_data, init_sky, plot_grid_map, plot_healpix_map, plot_sky_circles, plot_sky_binned
from matplotlib import pyplot as plt
import numpy as np

from generate_mask_footprint_shapely import load_boundary_multipolygon, get_boundary_coordinates

# Make empty mollweide projection plot, with equatorial coordinates
ax = init_sky(galactic_plane_color=None, ecliptic_plane_color=None)

poly = load_boundary_multipolygon()

for geom in poly.geoms:
    xs, ys = geom.exterior.xy
    print(len(xs), len(ys))
    xs_rad = ax.projection_ra(np.array(xs))
    ys_rad = ax.projection_dec(np.array(ys))
    ax.fill(xs_rad, ys_rad, alpha=0.5)

#plt.savefig("moll.pdf")
plt.show()
