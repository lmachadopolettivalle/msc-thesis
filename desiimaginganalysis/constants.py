BANDS = ["g", "r", "z"]

BASS_MzLS = "BASS-MzLS"
DECaLS_NGC = "DECaLS-NGC"
DECaLS_SGC = "DECaLS-SGC"

ALL_REGIONS = (BASS_MzLS, DECaLS_NGC, DECaLS_SGC)

def map_region_to_north_south(region=BASS_MzLS):
    if region == BASS_MzLS:
        return "north"
    if (region == DECaLS_NGC) or (region == DECaLS_SGC):
        return "south"

    raise ValueError("Invalid region!")

BGS_BRIGHT = 1
BGS_FAINT = 0
# Note some targets may belong to WISE. We will not use these here.
BGS_WISE = -1
