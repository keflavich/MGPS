import pyavm
from astropy import units as u
import aplpy
import numpy as np
import copy
import matplotlib
from astropy.io import fits
from astropy import wcs
import reproject
import pylab as pl
import PIL
from PIL import ImageEnhance

import os

#b1 = fits.open('other_data/AG-Laboca-Planck.49.5.fits')
b1 = fits.open('other_data/jps50_data.fits')[0]
b2 = fits.open('other_data/W51_90cm_CDBB.fits')[0]
b3 = fits.open('W51_feathered_img_PlanckCombined.fits')[0]

hdr = wcs.WCS(b3.header).celestial.to_header()
hdr['NAXIS'] = 2
hdr['NAXIS1'] = b3.header['NAXIS1']
hdr['NAXIS2'] = b3.header['NAXIS2']
hdr = fits.Header.fromtextfile('hdr4096.hdr')

b1d,_ = reproject.reproject_interp((b1.data.squeeze().T, wcs.WCS(b1.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)
assert np.isfinite(np.nanmax(b1d))
b2d,_ = reproject.reproject_interp((b2.data.squeeze(), wcs.WCS(b2.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)
b3d,_ = reproject.reproject_interp((b3.data.squeeze(), wcs.WCS(b3.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)


def linearize(x, xmin=None, xmax=None, truncate=True):
    if np.isscalar(x):
        return x
    else:
        if xmin is None:
            xmin = np.nanmin(x)
        if xmax is None:
            xmax = np.nanmax(x)
        if truncate:
            x = np.copy(x)
            x[x<xmin] = xmin
            x[x>xmax] = xmax
        return ((x-xmin)/(xmax-xmin))

def logscale(arr, logexp=3.0, toint=True, relinearize=True, **kwargs):
    linarr = linearize(arr, **kwargs)
    if logexp is None:
        logarr = linarr
    else:
        logarr = np.log10(linarr * 10**logexp + 1)
    if relinearize:
        return linearize(logarr)
    elif toint:
        lla = linearize(logarr)*255
        return lla.astype('uint8')
    else:
        return logarr

def expscale(arr, exp=2, toint=True, **kwargs):
    linarr = linearize(arr, **kwargs)
    if toint:
        lla = linearize(linarr**exp)*255
        return lla.astype('uint8')
    else:
        return linarr**exp


red = linearize(b2d, xmin=-0.05, xmax=0.4) * 0.75 + linearize(b2d, xmin=0.4, xmax=0.7)*0.25
#red = linearize(b2d, xmin=0.01, xmax=0.4)
#green = (linearize(b3d, 0, 4))
green = (logscale(b3d, xmin=0.003, xmax=4, logexp=4, toint=False))
green = (logscale(b3d, xmin=0.003, xmax=2, logexp=4, toint=False)) * 0.75 + linearize(b3d, xmin=2, xmax=8)*0.25
#blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
blue = logscale(b1d, xmin=-10, xmax=4000, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)
slices = slice(300,1300), slice(400,1600), slice(None) #[300:1300,400:1600]
slices = slice(None),slice(None),slice(None),
rgb_crop = rgb[slices]

mywcs = wcs.WCS(hdr)[slices[:2]]
avm = pyavm.AVM.from_wcs(mywcs)

outfn = "W51_RGB_90cm_MGPS_ATLASGAL.png"
im = PIL.Image.fromarray((rgb_crop*255).astype('uint8')[::-1,:])
im.save(outfn)
avm.embed(outfn, outfn)

pl.matplotlib.rc_file('pubfiguresrc')

pl.figure(2).clf()
FF = aplpy.FITSFigure(outfn, figure=pl.figure(2))
FF.show_rgb(outfn)
#FF.set_tick_labels_format('dd.d','dd.d')
FF.save('aplpy_'+outfn)

FF.add_scalebar(((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value)
FF.scalebar.set_label("10 pc")
FF.scalebar.set_font_size(18)
FF.scalebar.set_font_weight('bold')
FF.scalebar.set_color('w')
FF.scalebar.set_linewidth(3)
FF.save('aplpy_scalebar_'+outfn)
