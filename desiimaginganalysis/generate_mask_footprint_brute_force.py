from astropy.io import fits

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
        print(len(footprint))

def mask(ra, dec):
    global footprint
    def is_in_brick(brick, r, d):
        return (r >= brick["ralow"]) and (r <= brick["rahigh"]) and (d >= brick["declow"]) and (d <= brick["dechigh"])

    for area in footprint:
        if is_in_brick(area, ra, dec):
            return True
    return False
