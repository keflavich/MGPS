"""
Redo the smoothing treating zero-pixels (those that have no data) as NaN.  This
is particularly important for G34.
"""
import numpy as np
from astropy.io import fits
import glob
from astropy.convolution import convolve_fft, Gaussian2DKernel

fwhmfactor = np.sqrt(8*np.log(2))

for fn in glob.glob("/Volumes/external/mgps/Feb5_2019/*_final.fits"):

    fh = fits.open(fn)

    data = fh[0].data
    data[data==0] = np.nan

    smdata = convolve_fft(data, Gaussian2DKernel(4/fwhmfactor), allow_huge=True)

    bmsize = (8.1**2 + 4**2)**0.5

    fh[0].data = smdata
    fh[0].header['BMAJ'] = bmsize/3600.
    fh[0].header['BMIN'] = bmsize/3600.
    fh.writeto(fn.replace("_final","_final_smooth4"), overwrite=True)
