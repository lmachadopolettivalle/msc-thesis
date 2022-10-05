from astropy.io import fits
from shapely.geometry import MultiPolygon, Polygon, box, Point
from shapely.ops import unary_union

import pickle

# After constructing and saving the footprint Polygon,
# this function loads the savef footprint
# into a polygon object, and returns it
def load_boundary_multipolygon():
    ### Load pickled boundary
    with open("footprint.pickle", "rb") as f:
        poly = pickle.load(f)

    return poly

# After constructing and saving the footprint Polygon,
# this function can determine whether a given (RA, Dec)
# is within the footprint
def mask(ra, dec):
    boundary_polygon = load_boundary_multipolygon()

    point = Point(ra, dec)

    return point.within(boundary_polygon)

# If this file is invoked with python,
# generate and save the boundary polygon
# based on the DR9 survey brick files
if __name__ == "__main__":
    PATH_TO_FOOTPRINT_FILES = "/cluster/home/lmachado/msc-thesis/dr9/footprints/"
    FOOTPRINT_FILES = ["survey-bricks-dr9-north.fits", "survey-bricks-dr9-south.fits"]

    footprint = []

    for filename in FOOTPRINT_FILES:
        with fits.open(f"{PATH_TO_FOOTPRINT_FILES}/{filename}") as f:
            data = f[1].data

            footprint.extend([
                {
                    "ralow": i["ra1"],
                    "rahigh": i["ra2"],
                    "declow": i["dec1"],
                    "dechigh": i["dec2"],
                }
                for i in data
            ])

    L = [
        box(i["ralow"], i["declow"], i["rahigh"], i["dechigh"])
        for i in footprint
    ]


    P = unary_union(L)
    if P.geom_type == 'Polygon':
        P = MultiPolygon([P])

    ### Pickle the boundary
    with open("footprint.pickle", "wb") as f:
        pickle.dump(P, f, pickle.HIGHEST_PROTOCOL)
