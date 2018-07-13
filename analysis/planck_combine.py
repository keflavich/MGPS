import os
from uvcombine import feather_simple, feather_plot
from astroquery.skyview import SkyView
from astropy.io import fits
from astropy import units as u
from astropy import coordinates
from astropy import wcs

files = {'G29':'G29/GAL_029_precon_2_arcsec_galactic_white_alt_svd_pass_8.fits',
         'G31':'GAL_031/GAL_031_precon_2_arcsec_pass_9.fits',
         'G01':'SgrB2/SgrB2_precon_2_arcsec_pass_8.fits',
         'G12':'W33/W33_precon_2_arcsec_cached_pass_19.fits',
         'G43':'W49/W49_precon_2_arcsec_pass_9.fits',
         'G49':'W51/W51_feathered_img.fits',
        }
dirs = {'G29':'G29',
        'G31':'GAL_031',
        'G01':'SgrB2',
        'G12':'W33',
        'G43':'W49',
        'G49':'W51',
       }

for regname, fn in files.items():

    outfile = "{0}_PlanckCombined.fits".format(fn.replace(".fits",""))
    #if os.path.exists(outfile):
    #    continue

    print(regname, fn)

    fh = fits.open(fn)
    ww = wcs.WCS(fh[0].header)
    center = coordinates.SkyCoord(fh[0].header['CRVAL1'], fh[0].header['CRVAL2'],
                                  frame=wcs.utils.wcs_to_celestial_frame(ww),
                                  unit=(u.deg, u.deg))

    planck_image = SkyView.get_images(center, 'Planck 100')[0][0]
    # bandpass information here: https://wiki.cosmos.esa.int/planckpla2015/index.php/The_RIMO#HFI_2
    # for now, just guess...
    planck_image.header['REFFREQ'] = 100e9
    planck_image.header['BMAJ'] = 9.65/60
    planck_image.header['BMIN'] = 9.65/60
    planck_image.header['BUNIT'] = 'K'
    planck_image.writeto('planck/{0}_100GHz.fits'.format(regname), overwrite=True)

    fh[0].header['BMAJ'] = 10/3600
    fh[0].header['BMIN'] = 10/3600
    fh[0].header['REFFREQ'] = 91.5e9
    fh[0].header['BUNIT'] = 'Jy/beam'

    rslt = feather_simple(fh[0], planck_image)

    fh[0].data = rslt.real
    fh.writeto(outfile, overwrite=True)
