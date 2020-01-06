import os
import numpy as np
import pylab as pl

from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import coordinates
from astropy import wcs
import aplpy

from files import files
from make_sed_cutout_image import make_sed_plot
from paths import catalog_figure_path, catalog_path
from utils import makesed, makeuplims
from catalog_flux_limits import flux_limits

import matplotlib as mpl

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.right'] = False
mpl.rcParams['xtick.color'] = 'w'
mpl.rcParams['ytick.color'] = 'w'


szinch = (15,12)
fig = pl.figure(5, figsize=szinch)
pl.pause(0.1)
for ii in range(5):
    fig.set_size_inches(szinch[0], szinch[1])
    pl.pause(0.1)
    try:
        assert np.all(fig.get_size_inches() == np.array(szinch))
        break
    except AssertionError:
        continue


for regname,fn in files.items():
    print(f"region {regname} file {fn}")
    for threshold,min_npix in ((4, 100), ):#(6, 15), (8, 15), (10, 15)):
        for min_delta in (1, ): #2):

            mgps_fn = fn

            catfn = f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch_gaussfits.ipac'
            if not os.path.exists(catfn):
                # during debugging stages, this can happen...
                print(f"Skipped {catfn}")
                continue
            else:
                print(f"Working on {catfn}")
            ppcat = Table.read(catfn, format='ascii.ipac')

            fig5 = pl.figure(5, figsize=szinch)
            assert np.all(fig.get_size_inches() == np.array(szinch))

            fig2 = pl.figure(2)
            fig2.clf()
            mgpsdetected = ppcat['rejected'] == 0
            for ii,row in enumerate(ppcat[mgpsdetected]):
                name = f'{regname}_{row["SourceName"].strip()}'

                # Disable this if redos are desired
                #if os.path.exists(f'{catalog_figure_path}/seds/SED_plot_{name}.png'):
                #    continue
                #if ii/4+1 > 9:
                #    break
                #ax = fig2.add_subplot(3, 3, int(ii/4) + 1)
                x,y = makesed(row)
                _,uplims = makeuplims(row)
                #ax.loglog(x, y, 'o-')

                frame = wcs.utils.wcs_to_celestial_frame(wcs.WCS(fits.getheader(mgps_fn)))

                crd = coordinates.SkyCoord(*row['x_cen', 'y_cen'], frame=frame.name, unit=(u.deg, u.deg))
                make_sed_plot(crd, mgps_fn, figure=fig5, regname=regname)

                # HACK: if it already exists, the labels are broken
                ax = fig5.add_subplot(4, 5, 20)
                if ax in fig5.axes:
                    ax.remove()
                if ax in fig5.axes:
                    fig5.axes.remove(ax)
                ax = fig5.add_subplot(4, 5, 20)
                ax.cla()
                ax.loglog(x, y, 'o-')
                ax.loglog(x, uplims, marker='v', color='orange', linestyle='none')
                ax.set_aspect('equal', 'box')
                ax.set_xlabel("Wavelength ($\mu$m)")
                ax.set_ylabel("Flux Density [mJy]")
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ax.set_xlim(10, 1e6)
                ax.xaxis.set_ticks([100,1e4,1e6])
                ax.xaxis.set_ticklabels(['100 $\mu$m', '1 cm', '1 m'])

                ax.xaxis.set_tick_params(color='k', labelcolor='k')
                ax.yaxis.set_tick_params(color='k', labelcolor='k')

                fig5.savefig(f'{catalog_figure_path}/seds/SED_plot_{name}.pdf', bbox_inches='tight')
                print(f"finished {name}")

# DEBUG                break
# DEBUG            break
# DEBUG        break
# DEBUG    break
# DEBUG
# DEBUGif False:

            #galhdr = fits.Header.fromtextfile('../GAL_031/g31gal.hdr')
            #data,_ = reproject.reproject_interp(files['G31'], galhdr)
            #hdu = fits.PrimaryHDU(data=data, header=galhdr)

            try:
                pl.figure(4, figsize=(12,12)).clf()
                FF = aplpy.FITSFigure(fn, figure=pl.figure(4), convention='calabretta')
                FF.show_grayscale()
                FF.tick_labels.set_font(color='k')

                ppcat['x_cen'].unit = u.deg
                ppcat['y_cen'].unit = u.deg
                #coords = coordinates.SkyCoord(ppcat['x_cen'].quantity, ppcat['y_cen'].quantity, frame='icrs').galactic
                coords = coordinates.SkyCoord(ppcat['x_cen'].quantity, ppcat['y_cen'].quantity, frame='galactic')

                herscheldetected = ppcat['HerschelDetected'] == 'True'
                spitzerdetected = ppcat['SpitzerDetected'] == 'True'
                cmdetected = ppcat['cmDetected'] == 'True'
                mgpsdetected = ppcat['rejected'] == 0
                mask = ((cmdetected) & (mgpsdetected))
                if any(mask):
                    FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='w', facecolor='w', marker='x', linewidth=1)
                mask = ((herscheldetected) & (mgpsdetected))
                if any(mask):
                    FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='c', facecolor='c', marker='+', linewidth=1)
                mask = ((cmdetected) & (~herscheldetected) & (mgpsdetected))
                if any(mask):
                    FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='r')
                mask = ((~cmdetected) & (~herscheldetected) & (mgpsdetected))
                if any(mask):
                    FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='lime')

                FF.savefig(f'{catalog_figure_path}/{regname}_catalog_overlay.pdf')
            except Exception as ex:
                print(ex)
