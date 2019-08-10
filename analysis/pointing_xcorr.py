import numpy as np
import os
import requests
import regions
import json
import shutil
import image_registration
from scipy import ndimage
from astroquery.magpis import Magpis
#from astroquery.higal import HiGal
import astroquery.exceptions
from astropy.wcs import WCS, utils as wcsutils
from astropy.io import fits
from astropy import units as u, coordinates
from astropy import convolution
from astropy.utils.data import download_file
from astropy import stats
from astropy.visualization import (MinMaxInterval, AsinhStretch,
                                   PercentileInterval,
                                   ImageNormalize)

import radio_beam

import pylab as pl
import matplotlib

import warnings
warnings.filterwarnings('ignore')

from paths import diagnostic_figure_path, pilotpaperpath, catalog_path
from files import files

Magpis.cache_location = '/Volumes/external/mgps/cache/'
#HiGal.cache_location = '/Volumes/external/mgps/cache/'

print("WARNING: This script requires big downloads.")

offset = {}

bolocam_base_url = 'https://irsa.ipac.caltech.edu/data/BOLOCAM_GPS/images/v2/INNER_GALAXY/maps/'
bolocam_filename_map = {
    'G01': 'v2.1_ds2_l001_13pca_map20.fits',
    'G12': 'v2.1_ds2_l012_13pca_map20.fits',
    'G29': 'v2.1_ds2_l030_13pca_map20_crop.fits',
    'G31': 'v2.1_ds2_l031_13pca_map20_crop.fits',
    'G34': 'v2.1_ds2_l035_13pca_map20_crop.fits',
    'G43': 'v2.1_ds2_l040_13pca_map20.fits',
    'G49': 'v2.1_ds2_l050_13pca_map20_crop.fits',
}
# override central coordinate for MGPS cutout because sometimes they're too small and exclude the important flux
center_coordinate = {
    'G49': coordinates.SkyCoord(49.5, -0.4, frame='galactic', unit=(u.deg, u.deg)),
    'G12': coordinates.SkyCoord(12.7, -0.15, frame='galactic', unit=(u.deg, u.deg)),
}
gps20_override = {
    'G49': '/Users/adam/work/w51/vla_old/W51-LBAND_Carray.fits',
    #'G49': '/Users/adam/work/w51/vla_old/W51-LBAND-feathered_ABCD.fits',
    'G01': '/Users/adam/work/gc/20cm_0.fits',
}

def download_bolocam_file(region, dlpath=Magpis.cache_location):
    expected_fpath = os.path.join(dlpath, bolocam_filename_map[region])
    if os.path.exists(expected_fpath):
        return expected_fpath
    else:
        # doesn't do anything, downloads are always to /var/tmp
        #os.chdir(dlpath)

        url = f"{bolocam_base_url}/{bolocam_filename_map[region]}"
        print(f"Downloading URL {url}")
        fpath = download_file(url)
        print(f"Downloaded to {fpath} and moved to {expected_fpath}")

        shutil.move(fpath, expected_fpath)

        return expected_fpath


beams = {'bolocam': 33*u.arcsec,
         'atlasgal': 20*u.arcsec,}
mgps_beam = radio_beam.Beam(10*u.arcsec)

for regname,fn in files.items():
#DEBUG for regname,fn in (('G01',files['G01']),):
#DEBUG for regname,fn in (('G49',files['G49']),):
    fh = fits.open(fn)[0]

    # calculate offsets against original pointing before
    # "apply_pointing_corrections" is applied
    if 'CRVAL1A' in fh.header:
        fh.header['CRVAL1'] = fh.header['CRVAL1A']
        fh.header['CRVAL2'] = fh.header['CRVAL2A']

    if regname == 'G34':
        # fix bad stripe by removing it
        fh.data[:,2052:2060] = np.nan
        fh.data[2125:2140,2045:2060] = np.nan

    offset[regname] = {}

    ww = WCS(fh.header)
    center = ww.wcs_pix2world(fh.data.shape[1]/2, fh.data.shape[0]/2, 0)
    mgps_pixscale = wcsutils.proj_plane_pixel_area(ww)**0.5*u.deg

    if regname in center_coordinate:
        coordinate = center_coordinate[regname]
    else:
        coordinate = coordinates.SkyCoord(center[0], center[1], unit=(u.deg, u.deg),
                                          frame=wcsutils.wcs_to_celestial_frame(ww))
    radius = np.max(fh.data.shape) * wcsutils.proj_plane_pixel_area(ww)**0.5*u.deg
    if radius > 1.25*u.deg:
        radius = 1.25*u.deg

    print(f"region {regname} file {fn}")
    for survey in Magpis.list_surveys():
        #if not ('atlasgal' in survey or 'bolocam' in survey):
        if survey not in ('gps20new',):
            continue
        offset[regname][survey] = {}

        if 'bolocam' in survey:
            bolocam_fn = download_bolocam_file(regname)
            print(f"Loading {bolocam_fn}")
            data = fits.open(bolocam_fn)
        else:
            magpis_fn = os.path.join(Magpis.cache_location, f'{survey}_{regname}_MAGPIS.fits')
            if os.path.exists(magpis_fn):
                data = fits.open(magpis_fn)
                print(f"Loaded {magpis_fn} from disk")
            elif survey == 'gps20new' and regname in gps20_override:
                magpis_fn = gps20_override[regname]
                data = fits.open(magpis_fn)
                print(f"Loaded {magpis_fn} from disk")
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
        hdu.data = hdu.data.squeeze()
        hdu.header = magpis_ww.to_header()

        assert hdu.header['NAXIS'] == 2
        assert hdu.data.ndim == 2

        pixscale = wcsutils.proj_plane_pixel_area(magpis_ww)**0.5*u.deg

        # convolve MGPS to atlasgal/bolocam resolution
        if survey in beams:
            convbeam = radio_beam.Beam(beams[survey]).deconvolve(mgps_beam)
            print(f"{regname}: Convolving MGPS data to {survey} resoln with convolving beam {str(convbeam)}")
            convdata = convolution.convolve_fft(fh.data, convbeam.as_kernel(mgps_pixscale), allow_huge=True)
            convfh = fits.PrimaryHDU(data=convdata, header=fh.header)
        else:
            convfh = fh

        if survey == 'bolocam':
            print(f"{regname}: Projecting {survey} to MGPS")
            # project bolocam to MGPS, except with bigger pixels
            target_header = ww[::5,::5].to_header()
            target_header['NAXIS'] = 2
            target_header['NAXIS1'] = int(fh.header['NAXIS1'] / 5)
            target_header['NAXIS2'] = int(fh.header['NAXIS2'] / 5)
            pixscale = wcsutils.proj_plane_pixel_area(ww[::5,::5])**0.5*u.deg
            proj_image2, proj_image1, header = \
                    image_registration.FITS_tools.match_fits(convfh, hdu,
                                                             header=target_header,
                                                             return_header=True
                                                            )
            fits.PrimaryHDU(data=proj_image2, header=target_header).writeto(f"{Magpis.cache_location}/MGPS_{regname}_smBolo.fits", overwrite=True)
        else:
            #if regname == 'G01':
            #    # special case for 20cm data...
            #    # flag out Sgr A
            #    hdu.data[2150-1781:2204-1781,2656-1534:2721-1534] = np.nan
            #    print("Removed Sgr A*")
            reg = regions.read_ds9(f'/Users/adam/work/mgps/regions/freefreemask_{regname}.reg')
            ww = WCS(hdu.header)

            pixscale1, pixscale2 = (wcs.utils.proj_plane_pixel_scales(ww),
                                    wcs.utils.proj_plane_pixel_scales(wcs.WCS(convfh.header)))

            # project MGPS to retrieved b/c retrieved is always smaller in MAGPIS case
            if all(pixscale1 < pixscale2):
                print(f"{regname}: Projecting MGPS to {survey} pix={pixscale1}")
                proj_image1, proj_image2, header = \
                        image_registration.FITS_tools.match_fits(hdu, convfh,
                                                                 return_header=True
                                                                )
            else:
                print(f"{regname}: Projecting {survey} to MGPS pix={pixscale2}")
                proj_image2, proj_image1, header = \
                        image_registration.FITS_tools.match_fits(convfh, hdu,
                                                                 return_header=True
                                                                )
                ww = wcs.WCS(convfh.header)

            pixreg = [rr.to_pixel(ww) for rr in reg]
            rmasks = [rr.to_mask() for rr in pixreg]
            mask = np.zeros(proj_image1.shape, dtype='bool')
            for rm in rmasks:
                mask[rm.bbox.slices] += rm.data.astype('bool')

            #if regname == 'G01':
            #    assert np.all(np.isnan(hdu.data[2150:2204,2656:2721]))
            #    assert np.all(np.isnan(proj_image1[2150:2204,2656:2721]))
            # just to be EXTRA sure...
            ok = np.isfinite(proj_image1) & np.isfinite(proj_image2)
            ok &= (proj_image1 > 0) & (proj_image2 > 0)
            proj_image1[~ok] = np.nan
            proj_image2[~ok] = np.nan
            proj_image1[~mask] = np.nan
            proj_image2[~mask] = np.nan

        #raise ValueError()

        xcorr = image_registration.chi2_shift(proj_image1, proj_image2,
                                              zeromean=False,
                                              return_error=True,
                                              upsample_factor=100.)

        offset[regname][survey]['nomeansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

        print(f"{regname}: MAGPIS {survey} = {xcorr} = {xcorr*pixscale.to(u.arcsec)}")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy,cr1,cr2,shf = xcorr
        #raise ValueError("STOP HERE")
        #dx,dy,wdx,wdy,edx,edy,ewdx,ewdy = xcorr
        xcorr = image_registration.chi2_shift(proj_image1, proj_image2,
                                              zeromean=True, return_error=True,
                                              upsample_factor=100.)
        print(f"{regname}: zero-mean MAGPIS {survey} = {xcorr} = {xcorr*pixscale.to(u.arcsec)}")

        offset[regname][survey]['meansub'] = xcorr, xcorr*pixscale.to(u.arcsec)

        # annoyingly necessary diagnostics
        pl.figure(1).clf()
        pl.subplot(2,2,1).imshow(proj_image1, origin='lower', norm=matplotlib.colors.LogNorm())
        #pl.subplot(2,2,1).contour(proj_image1>stats.mad_std(proj_image1, ignore_nan=True)*2, colors=['k']*2, levels=[0.5])
        pl.subplot(2,2,2).imshow(proj_image2, origin='lower', norm=matplotlib.colors.LogNorm())
        #pl.subplot(2,2,2).contour(proj_image2>stats.mad_std(proj_image2, ignore_nan=True)*2, colors=['k']*2, levels=[0.5])
        pky,pkx = np.unravel_index(np.nanargmax(proj_image1), proj_image1.shape)
        npix = 50
        if pky < npix:
            pky = npix
        if pkx < npix:
            pkx = npix
        pl.subplot(2,2,3).imshow(proj_image1[pky-npix:pky+npix,pkx-npix:pkx+npix], origin='lower', norm=matplotlib.colors.LogNorm())
        pl.subplot(2,2,3).contour(proj_image2[pky-npix:pky+npix,pkx-npix:pkx+npix], linewidths=[0.1]*10, colors=['k']*10)
        pl.subplot(2,2,4).imshow(proj_image2[pky-npix:pky+npix,pkx-npix:pkx+npix], origin='lower', norm=matplotlib.colors.LogNorm())
        pl.subplot(2,2,4).contour(proj_image1[pky-npix:pky+npix,pkx-npix:pkx+npix], linewidths=[0.1]*10, colors=['k']*10)
        pl.savefig(f'{diagnostic_figure_path}/{regname}_{survey}_xcorr_diagnostics.png', bbox_inches='tight')

        pl.figure(2).clf()
        diffim = (proj_image1/np.nanpercentile(proj_image1, 99) -
                  proj_image2/np.nanpercentile(proj_image2, 99))
        slices = ndimage.find_objects(np.isfinite(diffim))[0]
        diffim = (proj_image1/np.nanpercentile(proj_image1[slices], 99) -
                  proj_image2/np.nanpercentile(proj_image2[slices], 99))
        asinhnorm = ImageNormalize(diffim[slices],
                                   interval=PercentileInterval(99.5),
                                   stretch=AsinhStretch())

        pl.subplot(1,1,1).imshow(diffim[slices],
                                 origin='lower', norm=asinhnorm)

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



#print("bolocam zeromean")
#for reg in offset:
#    print(f"{reg:5s}: {offset[reg]['bolocam']['meansub'][1][:2]}")

print("gps20new zeromean")
for reg in offset:
    print(f"{reg:5s}: {offset[reg]['gps20new']['meansub'][1][:2]}")
#print("bolocam no zeromean")
#for reg in offset:
#    print(f"{reg:5s}: {offset[reg]['bolocam']['nomeansub'][1]}")
#print("atlasgal")
#for reg in offset:
#    print(f"{reg:5s}: {offset[reg]['atlasgal']['meansub'][1]}")

for key in offset:
    offset[key]['bolocam']['meansub'] = offset[key]['bolocam']['meansub'][0],list(offset[key]['bolocam']['meansub'][1].value)
    offset[key]['bolocam']['nomeansub'] = offset[key]['bolocam']['nomeansub'][0],list(offset[key]['bolocam']['nomeansub'][1].value)
    offset[key]['gps20new']['meansub'] = offset[key]['gps20new']['meansub'][0],list(offset[key]['gps20new']['meansub'][1].value)
    offset[key]['gps20new']['nomeansub'] = offset[key]['gps20new']['nomeansub'][0],list(offset[key]['gps20new']['nomeansub'][1].value)

if len(offset) >= 6:
    # don't write (overwrite) the data unless all 6 fields have been fit
    # (allows running this script in "DEBUG" mode")
    with open(f"{catalog_path}/position_offsets.json", "w") as fh:
        json.dump(offset, fh)

#with open(f"{pilotpaperpath}/position_offsets_table.tex", "w") as fh:
#    fh.write(

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
