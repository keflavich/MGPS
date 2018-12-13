import numpy as np
import os
import glob
from astropy.io import fits
from astropy import stats
from astropy import wcs
import regions

import warnings
warnings.simplefilter('ignore', wcs.FITSFixedWarning)

dirs = [
"G29",
"SgrB2",
"W33",
"W43",
"W49",
"W51",
]

for dir in dirs:
    for fn in glob.glob("../{0}/*.fits".format(dir)):
        if 'ALMA' in fn:
            continue
        regionfn = os.path.split(fn)[0]+"/footprint.reg"
        region = regions.io.read_ds9(regionfn)[0]

        fh = fits.open(fn)
        header = fh[0].header
        ww = wcs.WCS(header)

        reg = region.to_pixel(ww)
        mask = reg.to_mask()
        data = mask.cutout(fh[0].data)
        data[~mask.data.astype('bool')] = np.nan

        print("{1:0.5g} {0}".format(fn, stats.mad_std(data, ignore_nan=True)))
