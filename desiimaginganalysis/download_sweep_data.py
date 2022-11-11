from bs4 import BeautifulSoup
from os.path import exists
import requests
import subprocess

BASEURL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/"

REGION = "north" # north or south

SWEEPFILESURL = f"{BASEURL}/{REGION}/sweep/9.0/"

OUTPUT_DIRECTORY = f"/home/ipa/refreg/data/DESI/desi_legacy_catalog/dr9/{REGION}/sweeps/"

# Fetch HTML with list of sweeps
# and extract list of .fits files
r = requests.get(SWEEPFILESURL)
soup = BeautifulSoup(r.content, "html.parser")
rows = soup.findAll("tr")
sweepfiles = []
for row in rows:
    try:
        fitsfile = row.findAll("a")[0].get("href")
        if fitsfile.endswith(".fits"):
            sweepfiles.append(fitsfile)
    except:
        continue

# Download .fits files
for sweepfile in sweepfiles:
    path = f"{SWEEPFILESURL}/{sweepfile}"
    # Don't re-download existing files
    if exists(f"{OUTPUT_DIRECTORY}/{sweepfile}"):
        print(f"{sweepfile} has already been downloaded, skipping...")
        continue

    subprocess.run(["wget", "--limit-rate=100m", f"--directory-prefix={OUTPUT_DIRECTORY}", path])
