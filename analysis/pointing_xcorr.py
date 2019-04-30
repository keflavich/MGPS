import numpy as np
import os
import requests
import image_registration
from astroquery.magpis import Magpis
from astroquery.higal import HiGal
import astroquery.exceptions
from astropy.wcs import WCS, utils as wcsutils
from astropy.io import fits
from astropy import units as u, coordinates

import warnings
warnings.filterwarnings('ignore')

from files import files

Magpis.cache_location = '/Volumes/external/mgps/cache/'
HiGal.cache_location = '/Volumes/external/mgps/cache/'

print("WARNING: This script requires big downloads.")
print("As of April 29, this didn't seem to work all that well.")

offset = {}

for regname,fn in files.items():
    fh = fits.open(fn)[0]

    offset[regname] = {}

    ww = WCS(fh.header)
    center = ww.wcs_pix2world(fh.data.shape[1]/2, fh.data.shape[0]/2, 0)

    coordinate = coordinates.SkyCoord(center[0], center[1], unit=(u.deg, u.deg),
                                      frame=wcsutils.wcs_to_celestial_frame(ww))
    radius = np.max(fh.data.shape) * wcsutils.proj_plane_pixel_area(ww)**0.5*u.deg
    if radius > 1*u.deg:
        radius = 1*u.deg

    print(f"region {regname} file {fn}")
    for survey in Magpis.list_surveys():
        offset[regname][survey] = {}
        magpis_fn = os.path.join(Magpis.cache_location, f'{survey}_{regname}_MAGPIS.fits')
        if os.path.exists(magpis_fn):
            data = fits.open(magpis_fn)
            #print(f"Loaded {magpis_fn} from disk")
        else:
            try:
                #print(f"Downloading {magpis_fn} from {survey} at {coordinate} +/- {radius}")
                data = Magpis.get_images(coordinate, image_size=radius,
                                         survey=survey)
            except astroquery.exceptions.InvalidQueryError:
                #print(f"Failed to retrieve {survey} {regname}")
                continue
            data.writeto(magpis_fn)

        hdu = data[0]
        # ditch everything but WCS to force celestial
        magpis_ww = WCS(hdu.header).celestial
        hdu.header = magpis_ww.to_header()

        pixscale = wcsutils.proj_plane_pixel_area(magpis_ww)**0.5*u.deg

        # project MGPS to retrieved b/c retrieved is always smaller
        proj_image1, proj_image2, header = \
                image_registration.FITS_tools.match_fits(hdu, fh,
                                                         return_header=True,
                                                         sigma_cut=2)

        #raise ValueError()

        xcorr = image_registration.chi2_shift(proj_image1, proj_image2, return_error=True, upsample_factor=10.)

        offset[regname][survey]['nomeansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

        print(f"MAGPIS {survey} = {xcorr}")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy,cr1,cr2,shf = xcorr
        #raise ValueError("STOP HERE")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy = xcorr
        xcorr = image_registration.chi2_shift(proj_image1, proj_image2, zeromean=True, upsample_factor=10.)
        print(f"zero-mean MAGPIS {survey} = {xcorr}")

        offset[regname][survey]['meansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

    try:
        HiGal.TIMEOUT = 300
        higal_ims = HiGal.get_images(coordinate, radius=radius)
        for hgim in higal_ims:
            xcorr = image_registration.FITS_tools.register_fits(fh, hgim)
            print(f"HiGal{hgim[0].header['WAVELEN']} = {xcorr}")
            xcorr = image_registration.FITS_tools.register_fits(fh, hgim, zeromean=True)
            print(f"zero-mean HiGal{hgim[0].header['WAVELEN']} = {xcorr}")
    except requests.exceptions.ReadTimeout:
        print("Skipped HiGal because of timeout.")
        continue
