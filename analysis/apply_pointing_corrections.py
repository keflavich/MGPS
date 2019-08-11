import json
from astropy import table
from astropy.io import fits
from astropy import log
from paths import catalog_path, pilotpaperpath
from files import files

with open(f"{catalog_path}/position_offsets.json", "r") as fh:
    offsets = json.load(fh)

for regname,fn in files.items():
    for suffix in ("", "_PlanckCombined"):
        sfn = fn.replace(".fits", suffix+".fits")
        fh = fits.open(sfn, mode='update')
        if 'CRVAL1A' not in fh[0].header:
            log.info(f"Correcting {fn} with {offsets[regname]['gps20new']['nomeansub'][1]}")
            fh[0].header['CRVAL1A'] = fh[0].header['CRVAL1']
            fh[0].header['CRVAL2A'] = fh[0].header['CRVAL2']
            # signs determined empirically.  I guess the sign swap is because of
            # the sign of pixscale?  i.e., cdelt is negative in RA
            fh[0].header['CRVAL1'] += offsets[regname]['gps20new']['nomeansub'][1][0] / 3600.
            fh[0].header['CRVAL2'] -= offsets[regname]['gps20new']['nomeansub'][1][1] / 3600.

            fh.flush()
        else:
            log.info(f"Skipping {fn}")
            log.info(f"CRVALnA = {fh[0].header['CRVAL1A']}, {fh[0].header['CRVAL2A']}")
            log.info(f"CRVALn  = {fh[0].header['CRVAL1']}, {fh[0].header['CRVAL2']}")
