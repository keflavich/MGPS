import numpy as np
import os
import pylab as pl
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

beam_map = {'atlasgal': 19.2*u.arcsec}


assumed_temperature = 25*u.K
assumed_dustbeta = 1.75


def make_hiidust_plot(coordinate, mgpsfile, width=1*u.arcmin,
                      surveys=['atlasgal'], figure=None,
                      regname='GAL_031'):

    mgps_fh = fits.open(mgpsfile)[0]
    frame = wcs.utils.wcs_to_celestial_frame(wcs.WCS(mgps_fh.header))

    coordname = "{0:06.3f}_{1:06.3f}".format(coordinate.galactic.l.deg,
                                             coordinate.galactic.b.deg)

    mgps_cutout = Cutout2D(mgps_fh.data, coordinate.transform_to(frame.name), size=width*2, wcs=wcs.WCS(mgps_fh.header))
    print(f"Retrieving MAGPIS data for {coordname} ({coordinate.to_string()} {coordinate.frame.name})")
    # we're treating 'width' as a radius elsewhere, here it's a full width
    images = {survey:getimg(coordinate, image_size=width*2, survey=survey) for survey in surveys}
    images = {x:y for x,y in images.items() if y is not None}
    #images['mgps'] = [mgps_cutout]

    # coordinate stuff so images can be reprojected to same frame
    ww = mgps_cutout.wcs.celestial
    mgps_pixscale = (wcs.utils.proj_plane_pixel_area(ww)*u.deg**2)**0.5

    if figure is None:
        figure = pl.gcf()
    figure.clf()

    for ii, (survey,img) in enumerate(images.items()):

        new_img = img[0].data
        if hasattr(img[0], 'header'):
            outwcs = wcs.WCS(img[0].header)
        else:
            outwcs = img[0].wcs

        tgt_bm = Beam(beam_map[survey])
        convbm = tgt_bm.deconvolve(mgps_beam)

        mgps_sm = convolution.convolve_fft(mgps_cutout.data, convbm.as_kernel(mgps_pixscale))
        mgps_reproj,_ = reproject.reproject_interp((mgps_sm, mgps_cutout.wcs), outwcs, shape_out=img[0].data.shape)

        dust_pred = dust_emissivity.blackbody.modified_blackbody(u.Quantity([wlmap[survey].to(u.GHz, u.spectral()), mustang_central_frequency]),
                                                                 assumed_temperature,
                                                                 beta=assumed_dustbeta)

        # assumes "surv" is dust
        surv_to_mgps = new_img * dust_pred[1]/dust_pred[0]
        print(f"agal to mgps ratio: {dust_pred[1]/dust_pred[0]}")

        dusty = surv_to_mgps / tgt_bm.sr.value
        freefree = (mgps_reproj / mgps_beam.sr.value - dusty)
        print(img[0].data.max(), mgps_sm.max())
        print(np.nanmax(dusty), np.nanmax(freefree), np.nanmax(mgps_reproj / mgps_beam.sr.value))

        norm = visualization.ImageNormalize(dusty,
                                            interval=visualization.ManualInterval(np.percentile(dusty, 0.05),
                                                                                  np.percentile(dusty, 99.95)),
                                            stretch=visualization.AsinhStretch(),
                                           )
        mgpsnorm = visualization.ImageNormalize(mgps_cutout.data,
                                                interval=visualization.PercentileInterval(99.95),
                                                stretch=visualization.AsinhStretch(),)

        ax0 = figure.add_subplot(1, 3, 1, projection=mgps_cutout.wcs)
        ax0.imshow(mgps_cutout.data, origin='lower', interpolation='none', norm=mgpsnorm)
        ax0.set_title("3 mm")
        ax1 = figure.add_subplot(1, 3, 2, projection=outwcs)
        ax1.imshow(dusty, origin='lower', interpolation='none', norm=norm)
        ax1.set_title("Dust")
        ax2 = figure.add_subplot(1, 3, 3, projection=outwcs)
        ax2.imshow(freefree, origin='lower', interpolation='none', norm=norm)
        ax2.set_title("Free-Free")

        for ax in (ax0, ax1, ax2):
            ax.set_xlabel("RA")
            ax.set_ylabel("Dec")
            ax.tick_params(direction='in')
            ax.tick_params(color='w')


        ax1.coords[1].set_axislabel("")
        ax1.coords[1].set_ticklabel_visible(False)
        ax2.coords[1].set_axislabel("")
        ax2.coords[1].set_ticklabel_visible(False)

        pl.subplots_adjust(hspace=0, wspace=0)

    pl.tight_layout()


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

    for reg in regs:

        for regname, mgpsfile in files.files.items():

            hdr = fits.getheader(mgpsfile)

            if hdr_contains_coord(hdr, reg.center):

                make_hiidust_plot(reg.center, mgpsfile, width=reg.radius, regname=regname)
                tgtname = reg.meta['label']

                pl.savefig(f'{paths.extended_figure_path}/{regname}_{tgtname}.pdf', bbox_inches='tight')
