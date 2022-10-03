from select_imaging_targets import *

from astropy.io import fits

target_count = 0

for filename in FILES:
    with fits.open(filename) as f:
        target_count += len(f[1].data)
        print(target_count)

print("total targets in sweep files: ", target_count)

print("selected targets: ", N)

print("fraction of objects that are kept by desitarget: ", N/target_count)
