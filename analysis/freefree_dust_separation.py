import numpy as np
import os
import reproject

from astroquery.magpis import Magpis
from astropy import units as u, coordinates
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
from astropy import visualization
from astropy import convolution
from radio_beam import Beam

# https://github.com/keflavich/dust_emissivity
import dust_emissivity

import paths
import files

from make_sed_cutout_image import wlmap, getimg
from constants import mgps_beam, mustang_central_frequency

import pylab as pl

import warnings
warnings.filterwarnings('ignore')

beam_map = {'atlasgal': 19.2*u.arcsec}


assumed_temperature = 25*u.K
assumed_dustbeta = 1.5


def make_hiidust_plot(reg, mgpsfile, width=1*u.arcmin,
                      surveys=['atlasgal'], figure=None,
                      regname='GAL_031',
                      fifth_panel_synchro=False,
                     ):

    mgps_fh = fits.open(mgpsfile)[0]
    frame = wcs.utils.wcs_to_celestial_frame(wcs.WCS(mgps_fh.header))

    coordinate = reg.center
    coordname = "{0:06.3f}_{1:06.3f}".format(coordinate.galactic.l.deg,
                                             coordinate.galactic.b.deg)

    mgps_cutout = Cutout2D(mgps_fh.data, coordinate.transform_to(frame.name), size=width*2, wcs=wcs.WCS(mgps_fh.header))
    print(f"Retrieving MAGPIS data for {coordname} ({coordinate.to_string()} {coordinate.frame.name})")
    # we're treating 'width' as a radius elsewhere, here it's a full width
    images = {survey:getimg(coordinate, image_size=width*2, survey=survey) for survey in surveys}
    images = {x:y for x,y in images.items() if y is not None}
    assert len(images) > 0
    #images['mgps'] = [mgps_cutout]

    # coordinate stuff so images can be reprojected to same frame
    ww = mgps_cutout.wcs.celestial
    mgps_pixscale = (wcs.utils.proj_plane_pixel_area(ww)*u.deg**2)**0.5

    if figure is None:
        figure = pl.gcf()
    figure.clf()

    (survey,img), = images.items()

    new_img = img[0].data
    if hasattr(img[0], 'header'):
        outwcs = wcs.WCS(img[0].header)
    else:
        outwcs = img[0].wcs

    agal_bm = tgt_bm = Beam(beam_map[survey])
    convbm = tgt_bm.deconvolve(mgps_beam)

    mgps_sm = convolution.convolve_fft(mgps_cutout.data, convbm.as_kernel(mgps_pixscale))
    mgps_reproj,_ = reproject.reproject_interp((mgps_sm, mgps_cutout.wcs), outwcs, shape_out=img[0].data.shape)

    mgpsjysr = mgps_cutout.data / mgps_beam.sr.value

    dust_pred = dust_emissivity.blackbody.modified_blackbody(u.Quantity([wlmap[survey].to(u.GHz,
                                                                                          u.spectral()),
                                                                         mustang_central_frequency]),
                                                             assumed_temperature,
                                                             beta=assumed_dustbeta)

    # assumes "surv" is dust
    surv_to_mgps = new_img * dust_pred[1]/dust_pred[0]
    print(f"{regname} {survey}")
    print(f"{survey} to mgps ratio: {dust_pred[1]/dust_pred[0]}")

    dusty = surv_to_mgps / tgt_bm.sr.value
    freefree = (mgps_reproj / mgps_beam.sr.value - dusty)
    print(img[0].data.max(), mgps_sm.max())
    print(np.nanmax(dusty), np.nanmax(freefree), np.nanmax(mgps_reproj / mgps_beam.sr.value))

    norm = visualization.ImageNormalize(freefree,
                                        interval=visualization.ManualInterval(np.nanpercentile(freefree, 0.1),
                                                                              np.nanpercentile(freefree, 99.9)),
                                        stretch=visualization.LogStretch(),
                                       )
    mgpsnorm = visualization.ImageNormalize(mgps_cutout.data,
                                            interval=visualization.PercentileInterval(99.95),
                                            stretch=visualization.LogStretch(),)
    print(f"interval: {norm.interval.vmin}, {norm.interval.vmax}")

    Magpis.cache_location = '/Volumes/external/mgps/cache/'

    ax0 = figure.add_subplot(1, 5, 3, projection=mgps_cutout.wcs)
    ax0.imshow(mgpsjysr, origin='lower', interpolation='none', norm=norm)
    ax0.set_title("3 mm")
    ax1 = figure.add_subplot(1, 5, 1, projection=outwcs)
    ax1.imshow(dusty, origin='lower', interpolation='none', norm=norm)
    ax1.set_title("870 $\\mu$m scaled")
    ax1.set_ylabel("Galactic Latitude")
    ax2 = figure.add_subplot(1, 5, 2, projection=outwcs)
    ax2.imshow(freefree, origin='lower', interpolation='none', norm=norm)
    ax2.set_title("3 mm Free-Free")

    for ax in (ax0, ax1, ax2):
        #ax.set_xlabel("Galactic Longitude")
        ax.tick_params(direction='in')
        ax.tick_params(color='w')


    ax0.coords[1].set_axislabel("")
    ax0.coords[1].set_ticklabel_visible(False)
    ax2.coords[1].set_axislabel("")
    ax2.coords[1].set_ticklabel_visible(False)

    pl.subplots_adjust(hspace=0, wspace=0)

    if 'G01' in regname:
        gps20im = fits.open('/Users/adam/work/gc/20cm_0.fits',)
    elif 'G49' in regname:
        gps20im = fits.open('/Users/adam/work/w51/vla_old/W51-LBAND-feathered_ABCD.fits')
        #gps20im = fits.open('/Users/adam/work/w51/vla_old/W51-LBAND_Carray.fits')
    else:
        gps20im = getimg(coordinate, image_size=width*2, survey='gps20new')

    reproj_gps20,_ = reproject.reproject_interp((gps20im[0].data.squeeze(),
                                                 wcs.WCS(gps20im[0].header).celestial),
                                                #mgps_fh.header)
    # refactoring to make a smaller cutout would make this faster....
                                                mgps_cutout.wcs,
                                                shape_out=mgps_cutout.data.shape)

    gps20cutout = Cutout2D(reproj_gps20, #gps20im[0].data.squeeze(),
                           coordinate.transform_to(frame.name), size=width*2,
                           wcs=mgps_cutout.wcs)
                           #wcs=wcs.WCS(mgps_fh.header))
                           #wcs.WCS(gps20im[0].header).celestial)
    ax3 = figure.add_subplot(1, 5, 5, projection=gps20cutout.wcs)

    gps20_bm = Beam.from_fits_header(gps20im[0].header)
    print(f"GPS 20 beam: {gps20_bm.__repr__()}")

    norm20 = visualization.ImageNormalize(gps20cutout.data,
                                        interval=visualization.ManualInterval(np.nanpercentile(gps20cutout.data, 0.5),
                                                                              np.nanpercentile(gps20cutout.data, 99.9)),
                                        stretch=visualization.LogStretch(),
                                       )

    # use 0.12 per Loren's suggestion
    freefree_20cm_to_3mm = (90*u.GHz/(1.4*u.GHz))**-0.12

    gps20_jysr = gps20cutout.data / gps20_bm.sr.value

    ax3.imshow(gps20_jysr * freefree_20cm_to_3mm, origin='lower', interpolation='none', norm=norm)
    ax3.set_title("20 cm scaled")
    ax0.set_xlabel("Galactic Longitude")
    ax3.coords[1].set_axislabel("")
    ax3.coords[1].set_ticklabel_visible(False)
    ax3.tick_params(direction='in')
    ax3.tick_params(color='w')



    # Fifth Panel:

    if fifth_panel_synchro:
        # MAGPIS data are high-resolution (comparable to but better than MGPS)
        # Zadeh data are low-resolution, 30ish arcsec
        # units: Jy/sr
        freefree_proj,_ = reproject.reproject_interp((freefree, outwcs),
                                                     gps20cutout.wcs,
                                                     shape_out=gps20cutout.data.shape)


        gps20_pixscale = (wcs.utils.proj_plane_pixel_area(gps20cutout.wcs)*u.deg**2)**0.5

        # depending on which image has higher resolution, convolve one to the other
        try:
            gps20convbm = tgt_bm.deconvolve(gps20_bm)
            gps20_jysr_sm = convolution.convolve_fft(gps20_jysr, gps20convbm.as_kernel(gps20_pixscale))
        except ValueError:
            gps20_jysr_sm = gps20_jysr
            ff_convbm = gps20_bm.deconvolve(tgt_bm)
            freefree_proj = convolution.convolve_fft(freefree_proj, ff_convbm.as_kernel(gps20_pixscale))

        ax4 = figure.add_subplot(1, 5, 5, projection=gps20cutout.wcs)

        # use the central frequency corresponding to an approximately flat spectrum (flat -> 89.72)
        freefree_3mm_to_20cm = 1/(90*u.GHz/(1.4*u.GHz))**-0.12
        #empirical_factor = 3 # freefree was coming out way too high, don't understand why yet
        synchro = gps20_jysr_sm - freefree_proj * freefree_3mm_to_20cm
        synchro[np.isnan(gps20_jysr) | (gps20_jysr == 0)] = np.nan

        synchroish_ratio = gps20_jysr_sm / (freefree_proj * freefree_3mm_to_20cm)

        #synchro = synchroish_ratio

        normsynchro = visualization.ImageNormalize(gps20_jysr_sm,
                                                   interval=visualization.ManualInterval(np.nanpercentile(gps20_jysr_sm,
                                                                                                          0.5),
                                                                                         np.nanpercentile(gps20_jysr_sm,
                                                                                                          99.9)),
                                                   stretch=visualization.LogStretch(),)

        ax4.imshow(synchro, origin='lower', interpolation='none', norm=normsynchro)
        ax4.set_title("Synchrotron")
        ax4.tick_params(direction='in')
        ax4.tick_params(color='w')
        ax4.coords[1].set_axislabel("")
        ax4.coords[1].set_ticklabel_visible(False)

        pl.tight_layout()
    else:
        # scale 20cm to match MGPS and subtract it

        gps20_pixscale = (wcs.utils.proj_plane_pixel_area(gps20cutout.wcs)*u.deg**2)**0.5


        if gps20_bm.sr < mgps_beam.sr:
            # smooth GPS20 to MGPS
            gps20convbm = mgps_beam.deconvolve(gps20_bm)
            gps20_jysr_sm = convolution.convolve_fft(gps20_jysr, gps20convbm.as_kernel(gps20_pixscale))
            gps20_jysr_sm[~np.isfinite(gps20_jysr)] = np.nan
            gps20_proj = gps20_jysr_sm
            #gps20_proj,_ = reproject.reproject_interp((gps20_jysr_sm, gps20cutout.wcs),
            #                                          ww,
            #                                          shape_out=mgps_cutout.data.shape)
        else:
            gps20_proj = gps20_jysr
            gps20_convbm = gps20_bm.deconvolve(mgps_beam)
            mgpsjysr = convolution.convolve_fft(mgpsjysr,
                                                gps20_convbm.as_kernel(mgps_pixscale))

        ax4 = figure.add_subplot(1, 5, 4, projection=mgps_cutout.wcs)

        # use the central frequency corresponding to an approximately flat spectrum (flat -> 89.72)
        dust20 = mgpsjysr - gps20_proj * freefree_20cm_to_3mm
        dust20[np.isnan(gps20_proj) | (gps20_proj == 0)] = np.nan

        normdust20 = visualization.ImageNormalize(mgpsjysr,
                                                  interval=visualization.ManualInterval(np.nanpercentile(mgpsjysr,
                                                                                                         0.5),
                                                                                        np.nanpercentile(mgpsjysr,
                                                                                                         99.9)),
                                                  stretch=visualization.LogStretch(),)

        # show smoothed 20 cm
        ax3.imshow(gps20_proj * freefree_20cm_to_3mm, origin='lower', interpolation='none', norm=norm)
        ax4.imshow(dust20, origin='lower', interpolation='none', norm=norm)
        ax4.set_title("3 mm Dust")
        ax4.tick_params(direction='in')
        ax4.tick_params(color='w')
        ax4.coords[1].set_axislabel("")
        ax4.coords[1].set_ticklabel_visible(False)

        pl.tight_layout()

    #elif 'G01' not in regname:
    #    norm.vmin = np.min([np.nanpercentile(dust20, 0.5), np.nanpercentile(freefree, 0.1)])
    if np.abs(np.nanpercentile(dust20, 0.5) - np.nanpercentile(freefree, 0.1)) < 1e8:
        norm.vmin = np.min([np.nanpercentile(dust20, 0.5), np.nanpercentile(freefree, 0.1)])
    if 'w49b' in reg.meta['text']:
        norm.vmin = np.min([np.nanpercentile(dust20, 8), np.nanpercentile(freefree, 0.1)])

    ax0.imshow(mgps_cutout.data / mgps_beam.sr.value, origin='lower', interpolation='none', norm=norm)
    ax1.imshow(dusty, origin='lower', interpolation='none', norm=norm)
    ax2.imshow(freefree, origin='lower', interpolation='none', norm=norm)
    ax3.imshow(gps20_proj * freefree_20cm_to_3mm, origin='lower', interpolation='none', norm=norm)
    ax4.imshow(dust20, origin='lower', interpolation='none', norm=norm)


if __name__ == "__main__":

    import regions

    regs = regions.read_ds9(os.path.join(paths.basepath, 'extended_regions.reg'))

    def hdr_contains_coord(header, coord):
        ww = wcs.WCS(header).celestial
        frame = wcs.utils.wcs_to_celestial_frame(ww)
        pixcrd = ww.wcs_world2pix(coord.transform_to(frame.name).spherical.lon,
                                  coord.transform_to(frame.name).spherical.lat,
                                  0
                                 )

        return pixcrd[0] > 0 and pixcrd[1] > 0 and pixcrd[0] < header['NAXIS1']-1 and pixcrd[1] < header['NAXIS2']-1

    pl.close(1)

    for reg in regs:

        for regname, mgpsfile in files.files.items():

            #mgpsfile = os.path.join('/Volumes/external/mgps/Feb5_2019/',
            mgpsfile = os.path.join('/Users/adam/work/mgps/Feb5_2019/',
                                    os.path.split(mgpsfile.replace(".fits","_PlanckCombined.fits"))[-1])

            ww = wcs.WCS(fits.getheader(mgpsfile))

            if ww.footprint_contains(reg.center):

                make_hiidust_plot(reg, mgpsfile, width=reg.radius, regname=regname,
                                  figure=pl.figure(1, figsize=(12,8)))
                tgtname = reg.meta['label']

                pl.savefig(f'{paths.extended_figure_path}/{regname}_{tgtname}_5panel.pdf', bbox_inches='tight')
