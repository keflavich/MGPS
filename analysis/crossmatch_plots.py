import os
import numpy as np
import pylab as pl

from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import coordinates
from astropy import wcs
import reproject
import aplpy

from files import files
from make_sed_cutout_image import make_sed_plot
from paths import catalog_figure_path, catalog_path

def makesed(row):
    frqscols = {#6*u.cm: 'Fint6cm_MAGPIS',
                6*u.cm: 'Fint6cm_CORNISH',
                20*u.cm: 'Fint20cm',
                3*u.mm: 'MUSTANG_dend_flux',
                350*u.um: 'Fint350um',
                250*u.um: 'Fint250um',
                160*u.um: 'Fint160um',
                70*u.um: 'Fint70um',
                24*u.um: 'Fint24um',
                870*u.um: 'Fint870um',
                1100*u.um: 'Fint1100um_40as',
                #8*u.um: 'Fint8um',
                #5.8*u.um: 'Fint5_8um',
                #4.5*u.um: 'Fint4_5um',
                #3.6*u.um: 'Fint3_6um',
               }
    freqs = sorted(frqscols.keys())
    values = [row.columns[frqscols[frq]].quantity[row.index] for frq in freqs]

    values = u.Quantity(values)
    values[values == 0] = np.nan
    return u.Quantity(freqs), values

flux_limits = {20*u.cm: 2*u.mJy, # MAGPIS: Helfand+ 2006
               6*u.cm: 2.5*u.mJy, # MAGPIS: Giveon+ 2005, CORNISH 2 mjy: Hoare+ 2012
              }


for regname,fn in files.items():
    for threshold,min_npix in ((4, 20), ):#(6, 15), (8, 15), (10, 15)):
        for min_delta in (1, ): #2):
            catfn = f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.ipac'
            if not os.path.exists(catfn):
                # during debugging stages, this can happen...
                continue
            ppcat = Table.read(catfn, format='ascii.ipac')

            herscheldetected = ppcat['HerschelDetected'] == 'True'
            spitzerdetected = ppcat['SpitzerDetected'] == 'True'
            cmdetected = ppcat['cmDetected'] == 'True'
            mgpsdetected = ppcat['rejected'] == 0

            print("Herschel-detected: {0} / {1} = {2}%".format((herscheldetected & mgpsdetected).sum(),
                                                               mgpsdetected.sum(),
                                                               100*(herscheldetected & mgpsdetected).sum() /
                                                               mgpsdetected.sum(),))
            print("Spitzer-detected: {0} / {1} = {2}%".format((spitzerdetected & mgpsdetected).sum(),
                                                               mgpsdetected.sum(),
                                                               100*(spitzerdetected & mgpsdetected).sum() /
                                                               mgpsdetected.sum(),))
            print("cm-detected: {0} / {1} = {2}%".format((cmdetected & mgpsdetected).sum(),
                                                               mgpsdetected.sum(),
                                                               100*(cmdetected & mgpsdetected).sum() /
                                                               mgpsdetected.sum(),))


            bins = np.logspace(-2,1,20)
            pl.figure(1)
            pl.clf()
            pl.subplot(3,1,1)
            pl.hist(ppcat['MUSTANG_dend_flux'][mgpsdetected], zorder=-5, histtype='step', color='k', bins=bins)
            pl.hist(ppcat['MUSTANG_dend_flux'][herscheldetected & mgpsdetected], bins=bins, label='Herschel-detected')
            pl.gca().set_xscale('log')
            pl.legend(loc='upper right', fontsize=12)
            pl.gca().set_xticklabels([])
            pl.subplot(3,1,2)
            pl.hist(ppcat['MUSTANG_dend_flux'][mgpsdetected], zorder=-5, histtype='step', color='k', bins=bins)
            pl.hist(ppcat['MUSTANG_dend_flux'][spitzerdetected & mgpsdetected], bins=bins, label='Spitzer-detected')
            pl.gca().set_xscale('log')
            pl.gca().set_xticklabels([])
            pl.legend(loc='upper right', fontsize=12)
            pl.subplot(3,1,3)
            pl.hist(ppcat['MUSTANG_dend_flux'][mgpsdetected], zorder=-5, histtype='step', color='k', bins=bins)
            pl.hist(ppcat['MUSTANG_dend_flux'][cmdetected & mgpsdetected], bins=bins, label='cm-detected')
            pl.xlabel("3 mm flux density [Jy]", fontsize=13)
            pl.gca().set_xscale('log')
            pl.legend(loc='upper right', fontsize=12)
            pl.subplots_adjust(hspace=0)

            pl.savefig(f'{catalog_figure_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_detection_histograms.pdf', bbox_inches='tight')

            hm = mgpsdetected & herscheldetected
            cm = mgpsdetected & cmdetected
            cmu = mgpsdetected & ~cmdetected

            ppcat['Fpeak6cm_MAGPIS'][mgpsdetected & (ppcat['Fpeak6cm_MAGPIS']==0)] = flux_limits[6*u.cm]
            ppcat['Fpeak20cm'][mgpsdetected & (ppcat['Fpeak20cm']==0)] = flux_limits[20*u.cm]

            pl.figure(3)
            pl.clf()
            pl.subplot(2,2,1)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm]/ppcat['Fint350um'][hm],
                      ppcat['Fint250um'][hm]/ppcat['Fint160um'][hm],
                      linestyle='none',
                      marker='o', alpha=0.75)
            pl.loglog(ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fint350um'][cm],
                      ppcat['Fint250um'][cm]/ppcat['Fint160um'][cm],
                      linestyle='none',
                      marker='o', alpha=0.5)
            pl.xlabel("3 mm / 350 $\mu$m")
            pl.ylabel("250 $\mu$m / 160 $\mu$m")

            pl.subplot(2,2,2)
            pl.loglog(ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fpeak6cm_MAGPIS'][cm]*1000,
                      ppcat['Fpeak6cm_MAGPIS'][cm]/ppcat['Fpeak20cm'][cm],
                      linestyle='none',
                      marker='o')
            #pl.loglog(ppcat['MUSTANG_dend_flux'][cmu]/ppcat['Fpeak6cm_MAGPIS'][cmu]*1000,
            #          ppcat['Fpeak6cm_MAGPIS'][cmu]/ppcat['Fpeak20cm'][cmu],
            #          linestyle='none',
            #          marker='^')
            pl.xlabel("3 mm / 6 cm")
            pl.ylabel("6 cm / 20 cm")

            pl.subplot(2,2,3)
            cm = mgpsdetected & cmdetected
            pl.loglog(ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fpeak20cm'][cm]*1000,
                    ppcat['Fpeak6cm_MAGPIS'][cm]/ppcat['Fint350um'][cm],
                    linestyle='none',
                    marker='o')
            pl.loglog(ppcat['MUSTANG_dend_flux'][cmu]/ppcat['Fpeak20cm'][cmu]*1000,
                    ppcat['Fpeak6cm_MAGPIS'][cmu]/ppcat['Fint350um'][cmu],
                    linestyle='none',
                    marker='^')
            pl.xlabel("3 mm / 20 cm")
            pl.ylabel("6 cm / 350 $\mu$m")

            pl.subplot(2,2,4)
            cm = mgpsdetected & cmdetected
            pl.loglog(ppcat['Fpeak6cm_MAGPIS'][cm]/ppcat['Fpeak20cm'][cm],
                      ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fint350um'][cm],
                      linestyle='none',
                      marker='o')
            pl.xlabel("6 cm / 20 cm")
            pl.ylabel("3 mm / 350 $\mu$m")


            pl.subplots_adjust(hspace=0.3)

            break
        break
    break

fig5 = pl.figure(5, figsize=(15,12))

fig2 = pl.figure(2)
fig2.clf()
mgpsdetected = ppcat['rejected'] == 0
for ii,row in enumerate(ppcat[mgpsdetected]):
    #if ii/4+1 > 9:
    #    break
    #ax = fig2.add_subplot(3, 3, int(ii/4) + 1)
    x,y = makesed(row)
    #ax.loglog(x, y, 'o-')

    mgps_fn = '../GAL_031/GAL031_5pass_1_.0.2_10mJy_10mJy_final_smooth4.fits'
    frame = wcs.utils.wcs_to_celestial_frame(wcs.WCS(fits.getheader(mgps_fn)))

    crd = coordinates.SkyCoord(*row['x_cen', 'y_cen'], frame=frame.name, unit=(u.deg, u.deg))
    make_sed_plot(crd, mgps_fn, figure=fig5, regname='GAL_031')

    ax = fig5.add_subplot(4, 5, 20)
    ax.loglog(x, y, 'o-')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel("Wavelength ($\mu$m)")
    ax.set_ylabel("Flux Density [Jy]")
    name = 'G031_{0}'.format(row['_idx'])
    name = f'{row["SourceName"]}'
    fig5.savefig(f'{catalog_figure_path}/seds/SED_plot_{name}.png', bbox_inches='tight')
    print(f"finished {name}")


galhdr = fits.Header.fromtextfile('../GAL_031/g31gal.hdr')
data,_ = reproject.reproject_interp(files['G31'], galhdr)
hdu = fits.PrimaryHDU(data=data, header=galhdr)

pl.figure(4).clf()
FF = aplpy.FITSFigure(hdu, figure=pl.figure(4))
FF.show_grayscale()

ppcat['x_cen'].unit = u.deg
ppcat['y_cen'].unit = u.deg
coords = coordinates.SkyCoord(ppcat['x_cen'].quantity, ppcat['y_cen'].quantity, frame='icrs').galactic

mask = ((cmdetected) & (mgpsdetected))
FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='w', facecolor='w', marker='x', linewidth=1)
mask = ((herscheldetected) & (mgpsdetected))
FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='c', facecolor='c', marker='+', linewidth=1)
mask = ((cmdetected) & (~herscheldetected) & (mgpsdetected))
FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='r')
mask = ((~cmdetected) & (~herscheldetected) & (mgpsdetected))
FF.show_markers(coords.l[mask], coords.b[mask], edgecolor='lime')

FF.savefig(f'{catalog_figure_path}/W43_catalog_overlay.pdf')
