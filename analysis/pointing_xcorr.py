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
from astropy import convolution
import radio_beam

import warnings
warnings.filterwarnings('ignore')

from files import files

Magpis.cache_location = '/Volumes/external/mgps/cache/'
HiGal.cache_location = '/Volumes/external/mgps/cache/'

print("WARNING: This script requires big downloads.")
print("As of April 29, this didn't seem to work all that well.")

offset = {}

beams = {'bolocam': 33*u.arcsec,
         'atlasgal': 20*u.arcsec,}
mgps_beam = radio_beam.Beam(10*u.arcsec)

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
        if not ('atlasgal' in survey or 'bolocam' in survey):
            continue
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

        # convolve MGPS to atlasgal/bolocam resolution
        convbeam = radio_beam.Beam(beams[survey]).deconvolve(mgps_beam)
        convdata = convolution.convolve_fft(fh.data, convbeam.as_kernel(pixscale), allow_huge=True)
        convfh = fits.PrimaryHDU(data=convdata, header=fh.header)

        # project MGPS to retrieved b/c retrieved is always smaller
        proj_image1, proj_image2, header = \
                image_registration.FITS_tools.match_fits(hdu, convfh,
                                                         return_header=True,
                                                         sigma_cut=2)

        #raise ValueError()

        xcorr = image_registration.chi2_shift(proj_image1, proj_image2, return_error=True, upsample_factor=100.)

        offset[regname][survey]['nomeansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

        print(f"MAGPIS {survey} = {xcorr}")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy,cr1,cr2,shf = xcorr
        #raise ValueError("STOP HERE")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy = xcorr
        xcorr = image_registration.chi2_shift(proj_image1, proj_image2, zeromean=True, upsample_factor=100.)
        print(f"zero-mean MAGPIS {survey} = {xcorr}")

        offset[regname][survey]['meansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

    #try:
    #    HiGal.TIMEOUT = 300
    #    higal_ims = HiGal.get_images(coordinate, radius=radius)
    #    for hgim in higal_ims:
    #        xcorr = image_registration.FITS_tools.register_fits(fh, hgim)
    #        print(f"HiGal{hgim[0].header['WAVELEN']} = {xcorr}")
    #        xcorr = image_registration.FITS_tools.register_fits(fh, hgim, zeromean=True)
    #        print(f"zero-mean HiGal{hgim[0].header['WAVELEN']} = {xcorr}")
    #except requests.exceptions.ReadTimeout:
    #    print("Skipped HiGal because of timeout.")
    #    continue



print("bolocam")
for reg in offset:
    print(f"{reg:5s}: {offset[reg]['bolocam']['meansub'][1]}")
print("atlasgal")
for reg in offset:
    print(f"{reg:5s}: {offset[reg]['atlasgal']['meansub'][1]}")

"""
G31  : [ 8.388 -3.348  4.068  3.636] arcsec
G12  : [  7.452 -10.188   1.116   1.152] arcsec
G43  : [-3.708 -9.972  0.828  0.972] arcsec
G49  : [ 2.052 12.852  1.044  1.152] arcsec
G01  : [11.196 -2.916  0.648  0.504] arcsec
G29  : [-4.644  3.924  5.076  4.428] arcsec
G34  : [1.548 2.34  1.368 1.548] arcsec


bolocam
G31  : [ 8.388 -3.42   4.104  3.672] arcsec
G12  : [  7.38  -10.116   1.116   1.188] arcsec
G43  : [ -3.852 -10.044   0.864   1.008] arcsec
G49  : [ 2.124 12.924  1.008  1.152] arcsec
G01  : [11.268 -2.988  0.648  0.504] arcsec
G29  : [-4.644  3.852  5.112  4.5  ] arcsec
G34  : [1.404 2.196 1.404 1.62 ] arcsec
atlasgal
G31  : [-0.03000006  0.39000078  1.98000396  1.83000366] arcsec
G12  : [-2.19000438 -6.1500123   0.63000126  0.6000012 ] arcsec
G43  : [0.51000102 1.71000342 0.42000084 0.51000102] arcsec
G49  : [ 7.77001554 32.01006402  0.54000108  0.6000012 ] arcsec
G01  : [-5.73001146 -2.73000546  0.24000048  0.24000048] arcsec
G29  : [-6.87001374  0.03000006  2.28000456  2.2500045 ] arcsec
G34  : [-0.57000114  5.2500105   0.6000012   0.63000126] arcsec
"""
