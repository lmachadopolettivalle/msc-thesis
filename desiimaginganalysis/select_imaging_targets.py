from desitarget.cuts import select_targets
from desitarget.targets import desi_mask, bgs_mask

FILES = ["sweep-040p025-050p030.fits", "sweep-260p025-270p030.fits", "sweep-010p000-020p005.fits", "sweep-020p025-030p030.fits", "sweep-240p025-250p030.fits"]
#FILES = ["sweep-020m050-030m045.fits"]

# Select BGS targets from sweep files
targets, infiles = select_targets(
    FILES,
    numproc=1,
    tcnames=["BGS"],
    survey="main",
    backup=False,
    return_infiles=True
)

# Load RA, DEC from BGS targets
RAs = [i["RA"] for i in targets]
DECs = [i["DEC"] for i in targets]
N = len(targets)

# Determine to which BGS a given target belongs
print(bgs_mask.names(targets[0]["BGS_TARGET"]))
