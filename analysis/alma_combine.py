import numpy as np
import uvcombine
import radio_beam
import reproject
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, visualization; from astropy.io import fits
import pylab as pl

almafn = 'SgrB2/SgrB2_selfcal_full_TCTE_selfcal5_ampphase_taylorterms_multiscale.image.tt0.pbcor.fits'
almafh = fits.open(almafn)[0]
loresfn = 'SgrB2/SgrB2_precon_2_arcsec_pass_9.fits'
loresfh_ = fits.open(loresfn)[0]

loresfwhm = 10*u.arcsec
loresbm = radio_beam.Beam(loresfwhm)

loresfh_header = loresfh_.header
loresfh_header.update(wcs.WCS(loresfh_.header)[1415:1815,958:1358].to_header())
loresfh = fits.PrimaryHDU(data=loresfh_.data[1415:1815,958:1358],
                          header=loresfh_header)

almabm = radio_beam.Beam.from_fits_header(almafh.header)

# hack for bad pixels?
#loresfh.data[1623:1629, 1157:1162] = np.nan
#loresfh.data[208:214, 199:204] = np.nan
pixscale_alma = wcs.utils.proj_plane_pixel_area(wcs.WCS(almafh.header))**0.5*u.deg
alma_in_gbt_units, almahdr_gbunits = uvcombine.uvcombine.match_flux_units(almafh.data, almafh.header, loresfh.header)
sm_alma = convolution.convolve_fft(alma_in_gbt_units, loresbm.deconvolve(almabm).as_kernel(pixscale_alma), allow_huge=True)
repr_alma,_ = reproject.reproject_interp((sm_alma, almafh.header), loresfh.header)
#loresfh.data[1623:1629, 1157:1162] = repr_alma[1623:1629, 1157:1162]
loresfh.data[208:214, 199:204] = repr_alma[208:214, 199:204]
loresfh.data[loresfh.data < repr_alma] = repr_alma[loresfh.data < repr_alma]



combined = uvcombine.feather_simple(almafn, lores=loresfh, lowresfwhm=loresfwhm, highresscalefactor=1)

hdr = fits.getheader(almafn)

fits.PrimaryHDU(data=combined.real,
                header=hdr).writeto('SgrB2/feathered_MGPS_ALMATCTE7m.fits',
                                    overwrite=True)

pl.figure(1)
rslts_thresh = uvcombine.feather_plot(almafn,
                                      lores='SgrB2/SgrB2_precon_2_arcsec_pass_9.fits',
                                      lowresfwhm=loresfwhm, hires_threshold=0.0005,
                                      lores_threshold=0.001)
pl.figure(2).clf()
pl.imshow(combined.real+0.01, origin='lower', interpolation='none', vmax=0.1, vmin=0.001,
          norm=pl.matplotlib.colors.LogNorm())
pl.axis((1620,2300,1842,2750))
