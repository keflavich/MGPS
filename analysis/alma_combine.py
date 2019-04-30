import numpy as np
import uvcombine
import radio_beam
import reproject
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, visualization; from astropy.io import fits
import pylab as pl

from constants import mustang_central_frequency, mustang_beam_fwhm

import paths

almafn = paths.root('SgrB2/SgrB2_selfcal_full_TCTE_selfcal5_ampphase_taylorterms_multiscale.image.tt0.pbcor.fits')
almafh = fits.open(almafn)[0]
loresfn = paths.root('SgrB2/SgrB2_precon_2_arcsec_pass_9.fits')
loresfn = paths.root('SgrB2/SgrB2_5pass_1_.0.2_10mJy_10mJy_final_smooth4.fits')
loresfn = '/Volumes/external/mgps/Feb5_2019/SgrB2_5pass_1_.0.2_10mJy_10mJy_w_session5_final_smooth4_PlanckCombined.fits'
loresfh_ = fits.open(loresfn)[0]

loresfwhm = mustang_beam_fwhm
loresbm = radio_beam.Beam(loresfwhm)

loresfh_header = loresfh_.header
if 'precon_2_arcsec_pass_9' in loresfn:
    loresfh_header.update(wcs.WCS(loresfh_.header)[1415:1815,958:1358].to_header())
    loresfh = fits.PrimaryHDU(data=loresfh_.data[1415:1815,958:1358],
                              header=loresfh_header)
elif '5pass_1_.0.2_10mJy' in loresfn:
    xc,yc = 1470,1990
    loresfh_header.update(wcs.WCS(loresfh_.header)[yc-350:yc+350,xc-350:xc+350].to_header())
    loresfh = fits.PrimaryHDU(data=loresfh_.data[yc-350:yc+350,xc-350:xc+350],
                              header=loresfh_header)

almabm = radio_beam.Beam.from_fits_header(almafh.header)

# hack for bad pixels?
#loresfh.data[1623:1629, 1157:1162] = np.nan
#loresfh.data[208:214, 199:204] = np.nan
pixscale_alma = wcs.utils.proj_plane_pixel_area(wcs.WCS(almafh.header))**0.5*u.deg
alma_in_gbt_units, almahdr_gbunits = uvcombine.uvcombine.match_flux_units(almafh.data, almafh.header, loresfh.header)
sm_alma = convolution.convolve_fft(alma_in_gbt_units, loresbm.deconvolve(almabm).as_kernel(pixscale_alma), allow_huge=True)
repr_alma,_ = reproject.reproject_interp((sm_alma, almafh.header), loresfh.header)
repr_gb,_ = reproject.reproject_interp(loresfh, almafh.header)
#loresfh.data[1623:1629, 1157:1162] = repr_alma[1623:1629, 1157:1162]
if 'precon_2_arcsec_pass_9' in loresfn:
    loresfh.data[208:214, 199:204] = repr_alma[208:214, 199:204]
    loresfh.data[loresfh.data < repr_alma] = repr_alma[loresfh.data < repr_alma]


alpha = 3
alma_freq = 96.34*u.GHz
mgps_freq = 92.44*u.GHz # alpha=3
lowresscalefactor = (alma_freq/mgps_freq)**alpha

# try to estimate the scalefactor empirically
pl.figure(4).clf()
pl.subplot(1,2,1)
sm_gb = convolution.convolve_fft(loresfh.data, convolution.Gaussian2DKernel(10))
unsharp_gb = loresfh.data - sm_gb
ok = (repr_alma > stats.mad_std(repr_alma, ignore_nan=True)*6) & (loresfh.data > stats.mad_std(loresfh.data, ignore_nan=True)*6)
pl.plot(repr_alma.flat, loresfh.data.flat, ',')
pl.plot(repr_alma.flat, unsharp_gb.flat, '.')
pl.plot(repr_alma[ok], unsharp_gb[ok], '.')

# TODO: do this using bayes_linear or ODR
scalefactor = np.nanmedian(unsharp_gb[ok] / repr_alma[ok])
pl.plot([0,7.5], [0,7.5 * scalefactor])
pl.subplot(2,2,2).imshow(loresfh.data, vmin=-0.1, vmax=1)
pl.subplot(2,2,4).imshow(loresfh.data - repr_alma * scalefactor, vmin=-0.1, vmax=1)


combined = uvcombine.feather_simple(almafn, lores=loresfh,
                                    lowresfwhm=loresfwhm, highresscalefactor=1,
                                    lowresscalefactor=1/scalefactor)

hdr = fits.getheader(almafn)

fits.PrimaryHDU(data=repr_alma, header=loresfh.header).writeto(paths.root('SgrB2/ALMAsmgridtoMGPS.fits'), overwrite=True)
fits.PrimaryHDU(data=sm_alma, header=almafh.header).writeto(paths.root('SgrB2/ALMAsmtoMGPS.fits'), overwrite=True)
loresfh.writeto(paths.root('SgrB2/MGPS_SgrB2_zoom_fixed.fits'), overwrite=True)


fits.PrimaryHDU(data=combined.real,
                header=hdr).writeto(paths.root('SgrB2/feathered_MGPS_ALMATCTE7m.fits'),
                                    overwrite=True)

pl.figure(1).clf()
rslts_thresh = uvcombine.feather_plot(almafn,
                                      lores=loresfn,
                                      lowresfwhm=loresfwhm, hires_threshold=0.0005,
                                      lores_threshold=0.001)
pl.figure(2).clf()
pl.imshow(combined.real+0.01, origin='lower', interpolation='none', vmax=0.1, vmin=0.001,
          norm=pl.matplotlib.colors.LogNorm())
pl.axis((1620,2300,1842,2750))


pl.figure(3).clf()
asinhnorm = lambda: visualization.ImageNormalize(visualization.AsinhStretch())

ax1 = pl.subplot(2,3,1)
im1 = ax1.imshow(combined.real+0.01, origin='lower', interpolation='none',
                 vmax=0.1, vmin=0.001, norm=asinhnorm())
ax1.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im1)

ax2 = pl.subplot(2,3,2)
im2 = ax2.imshow(almafh.data, origin='lower', interpolation='none', vmax=0.1,
                 vmin=-0.001, norm=asinhnorm())
ax2.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im2)

ax3 = pl.subplot(2,3,3)
im3 = ax3.imshow(combined.real - almafh.data, origin='lower',
                 interpolation='none', vmax=0.005, vmin=-0.005, norm=asinhnorm())
ax3.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im3)

ax4 = pl.subplot(2,3,4)
im4 = ax4.imshow(repr_gb, origin='lower', interpolation='none', vmax=4,
                 vmin=-0.005, norm=asinhnorm())
ax4.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im4)


ax5 = pl.subplot(2,3,5)
im5 = ax5.imshow(sm_alma * scalefactor, origin='lower',
                 interpolation='none', vmax=4, vmin=-0.005, norm=asinhnorm())
ax5.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im5)

ax6 = pl.subplot(2,3,6)
im6 = ax6.imshow(repr_gb - sm_alma * scalefactor, origin='lower', interpolation='none', vmax=0.1,
                 vmin=-0.1, norm=asinhnorm())
ax6.axis((1620,2300,1842,2750))
pl.colorbar(mappable=im6)
