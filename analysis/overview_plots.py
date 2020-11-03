import matplotlib
matplotlib.use('agg')
from dendrocat.aperture import Circle, Annulus
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy import units as u
from astropy import coordinates
from astropy.table import Column, Table
import regions
import pylab as pl
from paths import catalog_figure_path, catalog_path, overview_figure_path
from files import files

from constants import mustang_central_frequency, mustang_beam_fwhm

from astropy.visualization import (MinMaxInterval, AsinhStretch,
                                   PercentileInterval,
                                   ImageNormalize)


reglist = regions.io.read_ds9('cutout_regions.reg')
cutout_regions = {reg.meta['label']: reg for reg in reglist}

for regname,fn in files.items():
    for threshold,min_npix in ((4, 100),): # (6, 15), (8, 15), (10, 15)):
        for min_delta in (1, ):
            print(f"{regname}, {fn}")

            catalog = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch_gaussfits.ipac', format='ascii.ipac')

            fig = pl.figure(1)
            fig.clf()

            fh = fits.open(fn)
            data = fh[0].data
            header = fh[0].header
            # LONPOLE isn't very relevant and LATPOLE is not part of the coordinate
            # systems we're interested in.  From Calabretta 2002: "LATPOLEa is never
            # required for zenithal projections"
            try:
                del header['LONPOLE']
                del header['LATPOLE']
            except KeyError:
                pass
            ww = wcs.WCS(header)


            ax = fig.add_subplot(111, projection=ww)

            asinhnorm = ImageNormalize(data,
                                       interval=PercentileInterval(99.95),
                                       stretch=AsinhStretch())

            im = ax.imshow(data, cmap='gray_r', norm=asinhnorm)
            ax.set_xlabel("Galactic Longitude")
            ax.set_ylabel("Galactic Latitude")
            ax.tick_params(direction='in', color='k')


            #from mpl_toolkits.axes_grid1 import make_axes_locatable
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size="5%", pad=0.05)
            #cax.coords[0].grid(False)
            #cax.coords[1].grid(False)
            #cax.tick_params(direction='in')
            #cax.coords[0].set_ticks(alpha=0, color='w', size=0, values=[]*u.dimensionless_unscaled)
            #cax.coords[1].set_ticklabel_position('r')
            #cax.coords[1].set_axislabel_position('r')

            #pl.draw()

            #axbbox = ax.bbox.transformed(fig.transFigure.inverted())
            #print(f"Axis bbox in fig coords = {axbbox}")
            #axbbox = im.get_window_extent().transformed(fig.transFigure.inverted())
            #print(f"Image extent bbox in fig coords = {axbbox}")
            #cax = fig.add_axes([axbbox.x1 + 0.01, axbbox.y0, 0.03, axbbox.y1 - axbbox.y0])
            #assert cax is not None
            #cb = fig.colorbar(mappable=im, cax=cax)#, fraction=0.031, pad=0.04)
            #cb.set_label("$S_{3 mm}$ [Jy beam$^{-1}$]")
            #cax.set_ylabel("$S_{3 mm}$ [Jy beam$^{-1}$]")

            bbox = ax.get_position()
            bad_height = bbox.height
            print(f"bbox_height = {bbox.height}")
            fig.savefig(f"{overview_figure_path}/{regname}_overview.pdf", bbox_inches='tight')
            fig.savefig(f"{overview_figure_path}/{regname}_overview.png", bbox_inches='tight')
            bbox = ax.get_position()
            bad_height = bbox.height
            print(f"bbox_height = {bbox.height}")

            #ii = 0
            ## this is a painful hack to force the bbox to update
            #while bbox.height == bad_height:
            #    pl.pause(0.1)
            #    bbox = ax.get_position()
            #    print(f"bbox_height = {bbox.height}.  ii={ii}")
            #    ii += 1
            #    if ii > 10:
            #        break

            cax = fig.add_axes([bbox.x1+0.01, bbox.y0, 0.02, bbox.height])
            cb = fig.colorbar(mappable=im, cax=cax)
            cb.set_label("$S_{3 mm}$ [Jy beam$^{-1}$]")

            fig.savefig(f"{overview_figure_path}/{regname}_overview.pdf", bbox_inches='tight')
            fig.savefig(f"{overview_figure_path}/{regname}_overview.png", bbox_inches='tight')

            bolocamdetected = ~((catalog['Fint1100um'] == 0))
            atlasgaldetected = ~((catalog['Fint870um'] == 0))
            higaldetected = catalog['HerschelDetected'] == 'True'
            cmdetected = catalog['cmDetected'] == 'True'
            mmdetected = bolocamdetected | atlasgaldetected | higaldetected
            cm_mm_nondetection = (~cmdetected) & (~mmdetected)
            cm_no_mm_yes = (~cmdetected) & (mmdetected)



            kept = (catalog['rejected'] == 0) & (catalog['MorphologyClass'] == 'C')

            mask = kept & (mmdetected & cmdetected)
            if any(mask):
                keptpts, = ax.plot(catalog['x_cen'][mask], catalog['y_cen'][mask],
                                   marker='^', linestyle='none',
                                   markerfacecolor='none', markeredgecolor='m',
                                   transform=ax.get_transform('world'))
            mask = kept & (mmdetected & ~cmdetected)
            if any(mask):
                keptpts, = ax.plot(catalog['x_cen'][mask], catalog['y_cen'][mask],
                                   marker='s', linestyle='none',
                                   markerfacecolor='none', markeredgecolor='g',
                                   transform=ax.get_transform('world'))
            mask = kept & (cmdetected & ~mmdetected)
            if any(mask):
                keptpts, = ax.plot(catalog['x_cen'][mask], catalog['y_cen'][mask],
                                   marker='o', linestyle='none',
                                   markerfacecolor='none', markeredgecolor='b',
                                   transform=ax.get_transform('world'))
            mask = kept & cm_mm_nondetection
            if any(mask):
                keptpts, = ax.plot(catalog['x_cen'][mask], catalog['y_cen'][mask],
                                   marker='d', linestyle='none',
                                   markerfacecolor='none', markeredgecolor='orange',
                                   transform=ax.get_transform('world'))

            fig.savefig(f"{overview_figure_path}/{regname}_overview_withcatalog.pdf", bbox_inches='tight')
            fig.savefig(f"{overview_figure_path}/{regname}_overview_withcatalog.png", bbox_inches='tight')

            if any(~kept):
                allpts, = ax.plot(catalog['x_cen'][~kept], catalog['y_cen'][~kept],
                                  marker='x', linestyle='none',
                                  markerfacecolor='none', markeredgecolor='r',
                                  transform=ax.get_transform('world'))

            fig.savefig(f"{overview_figure_path}/{regname}_overview_withrejectcatalog.pdf", bbox_inches='tight')
            fig.savefig(f"{overview_figure_path}/{regname}_overview_withrejectcatalog.png", bbox_inches='tight')

            allpts.set_visible(False)

            hiicand_mask = catalog['HCHII_candidate'] == 'True'

            mask = kept & hiicand_mask
            if any(mask):
                hiicandpts, = ax.plot(catalog['x_cen'][mask],
                                      catalog['y_cen'][mask],
                                      marker='+', linestyle='none',
                                      markerfacecolor='b', markeredgecolor='c',
                                      transform=ax.get_transform('world'))


            fig.savefig(f"{overview_figure_path}/{regname}_overview_withcatalog_andHII.pdf", bbox_inches='tight')
            fig.savefig(f"{overview_figure_path}/{regname}_overview_withcatalog_andHII.png", bbox_inches='tight')

            axbbox = ax.bbox.transformed(fig.transFigure.inverted())
