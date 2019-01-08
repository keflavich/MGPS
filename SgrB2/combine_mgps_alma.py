"""
use analyis/alma_combine.py instead
"""
import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, visualization; from astropy.io import fits
import uvcombine
from uvcombine.uvcombine import regrid, file_in, match_flux_units
import pylab as pl

lores = mgps_image = 'SgrB2_5pass_1_.0.2_10mJy_10mJy_final_smooth4.fits'
hires = alma_image = 'SgrB2_selfcal_full_TCTE_selfcal5_ampphase_taylorterms_multiscale.image.tt0.pbcor.fits'

rslts_thresh = uvcombine.feather_plot(alma_image, lores=mgps_image,
                                      lowresfwhm=11*u.arcsec,
                                      hires_threshold=0.0005,
                                      lores_threshold=0.001)
pl.savefig('pspec_threshold.png')

rslts_nothresh = uvcombine.feather_plot(alma_image, lores=mgps_image, lowresfwhm=11*u.arcsec)
pl.savefig('pspec_nothreshold.png')

hdu_low, im_lowraw, header_low = file_in(lores)
hdu_hi, im_hi, header_hi = file_in(hires)




pl.figure(2).clf()
combined = uvcombine.feather_simple(alma_image, mgps_image, lowresfwhm=10*u.arcsec)

pl.imshow(visualization.AsinhStretch()(combined.real), vmax=0.1, origin='lower', interpolation='none')
hdr = fits.getheader(alma_image)

combined_default = uvcombine.feather_simple(alma_image, mgps_image)
fits.PrimaryHDU(data=np.abs(combined_default), header=hdr).writeto('feathered_MGPS_ALMATCTE7m.fits', overwrite=True)

combined_deconvsd = uvcombine.feather_simple(alma_image, mgps_image, deconvSD=True)
fits.PrimaryHDU(data=np.abs(combined_deconvsd), header=hdr).writeto('feathered_MGPS_ALMATCTE7m_deconvsd.fits', overwrite=True)

combined_lowres0p7 = uvcombine.feather_simple(alma_image, mgps_image, lowresscalefactor=0.7)
fits.PrimaryHDU(data=np.abs(combined_lowres0p7), header=hdr).writeto('feathered_MGPS_ALMATCTE7m_lowres0p7.fits', overwrite=True)

pl.figure(3).clf()
pl.subplot(2,2,1).imshow(visualization.AsinhStretch()(combined.real), vmax=0.1, origin='lower', interpolation='none')
pl.subplot(2,2,2).imshow(visualization.AsinhStretch()(combined_default.real), vmax=0.1, origin='lower', interpolation='none')
pl.subplot(2,2,3).imshow(visualization.AsinhStretch()(combined_deconvsd.real), vmax=0.1, origin='lower', interpolation='none')
pl.subplot(2,2,4).imshow(visualization.AsinhStretch()(combined_lowres0p7.real), vmax=0.1, origin='lower', interpolation='none')
