# File naming conventions are documented in https://www.legacysurvey.org/dr9/files/#tractor-catalogs-region-tractor
import subprocess

REGION = "south" # north or south

FILENAME = "sweep-040p025-050p030.fits"
PATH = f"https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/{REGION}/sweep/9.0/{FILENAME}"

subprocess.run(["wget", "--limit-rate=1m", PATH])
