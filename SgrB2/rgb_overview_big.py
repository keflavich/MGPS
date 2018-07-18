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
for b1,atlasgalscuba in ((fits.open('other_data/gc850_gal.fits')[0],'SCUBA'),
                         (fits.open('other_data/full_cmz_atlasgal.fits')[0], 'ATLASGAL')):
    #b2 = fits.open('other_data/VLA20cm_sgrb_12as.fits')[0]
    b2 = fits.open('other_data/20cm_0.fits')[0]
    b3p = fits.open('SgrB2_precon_2_arcsec_pass_8_PlanckCombined.fits')[0]
    b3 = fits.open('SgrB2_precon_2_arcsec_pass_8.fits')[0]

    hdr = fits.Header.fromtextfile('gal.hdr')

    b1d,_ = reproject.reproject_interp((b1.data.squeeze(), wcs.WCS(b1.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)
    assert np.isfinite(np.nanmax(b1d))
    b2d,_ = reproject.reproject_interp((b2.data.squeeze(), wcs.WCS(b2.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)
    b3d,_ = reproject.reproject_interp((b3.data.squeeze(), wcs.WCS(b3.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)
    b3dp,_ = reproject.reproject_interp((b3p.data.squeeze(), wcs.WCS(b3p.header).sub([wcs.WCSSUB_LONGITUDE, wcs.WCSSUB_LATITUDE])), hdr)


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


    red = linearize(b2d, xmin=0.025, xmax=0.5) * 0.75 + linearize(b2d, xmin=0.5, xmax=6)*0.25
    #red = linearize(b2d, xmin=0.01, xmax=0.4)
    #green = (linearize(b3d, 0, 4))
    #green = (logscale(b3d, xmin=0.003, xmax=4, logexp=4, toint=False))
    greenPlanck = (logscale(b3dp, xmin=0.02, xmax=0.2, logexp=2, toint=False)) * 0.75 + linearize(b3dp, xmin=0.2, xmax=3)*0.25
    green = (logscale(b3d, xmin=0.003, xmax=0.2, logexp=3, toint=False)) * 0.75 + linearize(b3d, xmin=0.2, xmax=3)*0.25
    #blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
    if atlasgalscuba == 'SCUBA':
        blue = (linearize(b1d, 0, 0.09) * 0.25 + logscale(b1d, xmin=0.09, xmax=1, logexp=4, toint=False)) * 0.75
        blue = linearize(b1d, xmin=0.06, xmax=0.3)
        blue = logscale(b1d, xmin=0.06, xmax=0.3, toint=False, logexp=1)
    else:
        blue = logscale(b1d, xmin=0.13, xmax=1, toint=False, logexp=1)*0.75 + linearize(b1d, xmin=1, xmax=15)*0.25
        blue = logscale(b1d, xmin=0.13, xmax=6, toint=False, logexp=1)*0.75 + linearize(b1d, xmin=6, xmax=30)*0.25

    rgb = np.array([red,green,blue]).T.swapaxes(0,1)
    slices = slice(None),slice(None),slice(None),
    slices = slice(400,-400),slice(None),slice(None),
    rgb_crop = rgb[slices]

    mywcs = wcs.WCS(hdr)[slices[:2]]
    avm = pyavm.AVM.from_wcs(mywcs)

    outfn = "SgrB2_RGB_20cm_MGPS_{0}.png".format(atlasgalscuba)
    im = PIL.Image.fromarray((rgb_crop*255).astype('uint8')[::-1,:])
    im.save(outfn)
    avm.embed(outfn, outfn)

    pl.matplotlib.rc_file('pubfiguresrc')

    pl.figure(2).clf()
    FF = aplpy.FITSFigure(outfn, figure=pl.figure(2))
    FF.show_rgb(outfn)
    #FF.set_tick_labels_format('dd.d','dd.d')
    FF.save('aplpy_'+outfn)

    FF.add_scalebar(((10*u.pc)/(6*u.kpc)*u.radian).to(u.deg).value)
    FF.scalebar.set_label("10 pc")
    FF.scalebar.set_font_size(18)
    FF.scalebar.set_font_weight('bold')
    FF.scalebar.set_color('w')
    FF.scalebar.set_linewidth(3)
    FF.save('aplpy_scalebar_'+outfn)




    rgb = np.array([red,greenPlanck,blue]).T.swapaxes(0,1)
    rgb_crop = rgb[slices]

    outfn = "SgrB2_RGB_20cm_MGPSplanck_{0}.png".format(atlasgalscuba)
    im = PIL.Image.fromarray((rgb_crop*255).astype('uint8')[::-1,:])
    im.save(outfn)
    avm.embed(outfn, outfn)

    pl.matplotlib.rc_file('pubfiguresrc')

    pl.figure(2).clf()
    FF = aplpy.FITSFigure(outfn, figure=pl.figure(2))
    FF.show_rgb(outfn)
    #FF.set_tick_labels_format('dd.d','dd.d')
    FF.save('aplpy_'+outfn)

    FF.add_scalebar(((10*u.pc)/(6*u.kpc)*u.radian).to(u.deg).value)
    FF.scalebar.set_label("10 pc")
    FF.scalebar.set_font_size(18)
    FF.scalebar.set_font_weight('bold')
    FF.scalebar.set_color('w')
    FF.scalebar.set_linewidth(3)
    FF.save('aplpy_scalebar_'+outfn)
