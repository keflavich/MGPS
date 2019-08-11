"""
Script to undo pointing corrections in case we need to re-do them
"""
from astropy.io import fits
from astropy import log
from files import files

for regname,fn in files.items():
    for suffix in ("", "_PlanckCombined"):
        sfn = fn.replace(".fits", suffix+".fits")
        fh = fits.open(sfn, mode='update')
        if 'CRVAL1A' in fh[0].header:
            log.info(f"unCorrecting {fn}")
            fh[0].header['CRVAL1'] = fh[0].header['CRVAL1A']
            fh[0].header['CRVAL2'] = fh[0].header['CRVAL2A']
            del fh[0].header['CRVAL1A']
            del fh[0].header['CRVAL2A']
            fh.flush()
        else:
            log.info(f"Skipping {fn}")
