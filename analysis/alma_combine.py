import uvcombine
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, visualization; from astropy.io import fits

almafn = 'SgrB2/SgrB2_selfcal_full_TCTE_selfcal5_ampphase_taylorterms_multiscale.image.tt0.pbcor.fits'

combined = uvcombine.feather_simple(almafn, lores='SgrB2/SgrB2_precon_2_arcsec_pass_9.fits', lowresfwhm=10*u.arcsec)

hdr = fits.getheader(almafn)

fits.PrimaryHDU(data=combined.real,
                header=hdr).writeto('SgrB2/feathered_MGPS_ALMATCTE7m.fits',
                                    overwrite=True)

rslts_thresh = uvcombine.feather_plot(almafn,
                                      lores='SgrB2/SgrB2_precon_2_arcsec_pass_9.fits',
                                      lowresfwhm=10*u.arcsec, hires_threshold=0.0005,
                                      lores_threshold=0.001)
