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
from utils import makesed

flux_limits = {20*u.cm: 2*u.mJy, # MAGPIS: Helfand+ 2006
               6*u.cm: 2.5*u.mJy, # MAGPIS: Giveon+ 2005, CORNISH 2 mjy: Hoare+ 2012
               70*u.um: 20*u.mJy, # from Fig 3 of Molinari 2016, ~20 MJy/sr rms w/6" beam -> 20 mJy/beam
               160*u.um: 26*u.mJy, # from Fig 3 of Molinari 2016, ~10 MJy/sr rms w/10" beam -> 26 mJy/beam
               350*u.um: 50*u.mJy, # total guess.  Close, though: 5 MJy/sr w/20" beam gives 53 mJy
               500*u.um: 85*u.mJy, # 2 MJy/sr in 4.24e-8 sr
               870*u.um: 70*u.mJy, # Csengeri+ 2014 fig 1
               1100*u.um: 50*u.mJy, # Ginsburg+ 2013, fig 1
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
            hmu = mgpsdetected & (~herscheldetected)
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

            # pl.subplot(2,2,2)
            # pl.loglog(ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fpeak6cm_MAGPIS'][cm]*1000,
            #           ppcat['Fpeak6cm_MAGPIS'][cm]/ppcat['Fpeak20cm'][cm],
            #           linestyle='none',
            #           marker='o')
            # #pl.loglog(ppcat['MUSTANG_dend_flux'][cmu]/ppcat['Fpeak6cm_MAGPIS'][cmu]*1000,
            # #          ppcat['Fpeak6cm_MAGPIS'][cmu]/ppcat['Fpeak20cm'][cmu],
            # #          linestyle='none',
            # #          marker='^')
            # pl.xlabel("3 mm / 6 cm")
            # pl.ylabel("6 cm / 20 cm")

            pl.subplot(2,2,2)
            # from http://herschel.esac.esa.int/hcss-doc-15.0/load/spire_drg/html/ch06s09.html
            f350um = (ppcat['Fpeak350um'].quantity * 1.95386e-8*u.sr).to(u.Jy)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm & cm]/ppcat['Fpeak6cm_MAGPIS'][hm & cm]*1000,
                      f350um[hm & cm]/ppcat['MUSTANG_dend_flux'][hm & cm],
                      linestyle='none',
                      marker='o', alpha=0.8)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hmu & cm]/ppcat['Fpeak20cm'][hmu & cm]*1000,
                      flux_limits[350*u.um].to(u.Jy)/ppcat['MUSTANG_dend_flux'][hmu & cm],
                      linestyle='none',
                      marker='v', alpha=0.8, zorder=-3)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm & cmu]/flux_limits[6*u.cm].to(u.Jy),
                      f350um[hm & cmu]/ppcat['MUSTANG_dend_flux'][hm & cmu],
                      linestyle='none',
                      marker='>', zorder=-5, alpha=0.8)
            xlims = pl.gca().get_xlim()
            ylims = pl.gca().get_ylim()
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**3]*2, 'k--', zorder=-10)
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**4]*2, 'k:', zorder=-10)
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**2]*2, 'k-', zorder=-10)
            #pl.plot([(6*u.cm / (3*u.mm)).decompose()**2]*2, [1,1e5], 'g--', zorder=-10)
            pl.plot([(6*u.cm / (3*u.mm)).decompose()**-0.1]*2, [1,1e5], 'r:', zorder=-10)
            pl.gca().set_xlim(*xlims)
            pl.gca().set_ylim(*ylims)
            pl.xlabel("3 mm / 6 cm")
            pl.ylabel("350 $\mu$m / 3 mm")


            pl.subplot(2,2,3)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm]/ppcat['Fpeak20cm'][hm]*1000,
                      ppcat['Fpeak70um'][hm]/ppcat['MUSTANG_dend_flux'][hm],
                      linestyle='none',
                      marker='o')
            pl.loglog(ppcat['MUSTANG_dend_flux'][hmu]/ppcat['Fpeak20cm'][hmu]*1000,
                      flux_limits[70*u.um]/ppcat['MUSTANG_dend_flux'][hmu],
                      linestyle='none',
                      marker='^')
            pl.xlabel("3 mm / 20 cm")
            pl.ylabel("70 $\mu$m / 3 mm")

            pl.subplot(2,2,4)
            cm = mgpsdetected & cmdetected
            pl.loglog(ppcat['Fpeak6cm_MAGPIS'][cm]/ppcat['Fpeak20cm'][cm],
                      ppcat['MUSTANG_dend_flux'][cm]/ppcat['Fint350um'][cm],
                      linestyle='none',
                      marker='o')
            pl.xlabel("6 cm / 20 cm")
            pl.ylabel("3 mm / 350 $\mu$m")


            pl.subplots_adjust(hspace=0.4, wspace=0.4)

            fign = f'{catalog_figure_path}/colorcolor_{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.pdf'
            pl.savefig(fign)


            pl.figure(4).clf()
            # from http://herschel.esac.esa.int/hcss-doc-15.0/load/spire_drg/html/ch06s09.html
            f350um = (ppcat['Fpeak350um'].quantity * 1.95386e-8*u.sr).to(u.Jy)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm & cm]/ppcat['Fpeak6cm_MAGPIS'][hm & cm]*1000,
                      f350um[hm & cm]/ppcat['MUSTANG_dend_flux'][hm & cm],
                      linestyle='none',
                      marker='o', alpha=0.8)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hmu & cm]/ppcat['Fpeak20cm'][hmu & cm]*1000,
                      flux_limits[350*u.um].to(u.Jy)/ppcat['MUSTANG_dend_flux'][hmu & cm],
                      linestyle='none',
                      marker='v', alpha=0.8, zorder=-3)

            # these are candidate HCHIIs: they have dust (because they're Hershel-detected)
            pl.loglog(ppcat['MUSTANG_dend_flux'][hm & cmu]/flux_limits[6*u.cm].to(u.Jy),
                      f350um[hm & cmu]/ppcat['MUSTANG_dend_flux'][hm & cmu],
                      linestyle='none',
                      marker='>', zorder=5, alpha=0.8)

            pl.loglog(ppcat['MUSTANG_dend_flux'][hmu & cmu]/flux_limits[6*u.cm].to(u.Jy),
                      flux_limits[350*u.um].to(u.Jy)/ppcat['MUSTANG_dend_flux'][hmu & cmu],
                      linestyle='none',
                      marker=[(0,0),(1,-1),(1,-0.25),(1,-1),(0.25,-1)], zorder=-5, alpha=0.8)

            xlims = pl.gca().get_xlim()
            ylims = pl.gca().get_ylim()
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**3]*2, 'k--', zorder=-10)
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**4]*2, 'k:', zorder=-10)
            pl.plot([1e-3, 1e3], [(3*u.mm/(350*u.um)).decompose().value**2]*2, 'k-', zorder=-10)
            #pl.plot([(6*u.cm / (3*u.mm)).decompose()**2]*2, [1,1e5], 'g--', zorder=-10)
            pl.plot([(6*u.cm / (3*u.mm)).decompose()**-0.1]*2, [1,1e5], 'r:', zorder=-10)
            pl.gca().set_xlim(*xlims)
            pl.gca().set_ylim(*ylims)
            pl.xlabel("3 mm / 6 cm")
            pl.ylabel("350 $\mu$m / 3 mm")

            fign = f'{catalog_figure_path}/colorcolor_3mm6cm_vs_350um3mm_{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.pdf'

            pl.figure(5).clf()
            # from http://herschel.esac.esa.int/hcss-doc-15.0/load/spire_drg/html/ch06s09.html
            f350um = (ppcat['Fpeak350um'].quantity * 1.95386e-8*u.sr).to(u.Jy)
            f500um = (ppcat['Fpeak500um'].quantity * 4.24e-8*u.sr).to(u.Jy)
            f870um = ppcat['Fint870um'].quantity.to(u.Jy) # (ppcat['Fpeak870um'].quantity * radio_beam.Beam(19.2*u.arcsec).sr).to(u.Jy)
            f1100um = ppcat['Fint1100um_40as'].quantity.to(u.Jy) * 1.46 # from Aguirre, Gins+ 2011

            pl.loglog(f350um[hm]/f1100um[hm],
                      f1100um[hm]/ppcat['MUSTANG_dend_flux'][hm].quantity.to(u.Jy),
                      linestyle='none',
                      marker='o', alpha=0.8)

            #pl.loglog([flux_limits[350*u.um].to(u.Jy) / flux_limits[1100*u.um].to(u.Jy)]*(hmu.sum()),
            #           flux_limits[1100*u.um].to(u.Jy)/ppcat['MUSTANG_dend_flux'][hmu].quantity.to(u.Jy),
            #           linestyle='none',
            #           marker=[(0,0),(1,-1),(1,-0.25),(1,-1),(0.25,-1)], zorder=-5, alpha=0.8)

            xlims = pl.gca().get_xlim()
            ylims = pl.gca().get_ylim()
            xlims = 9,101

            pl.plot([((350*u.um) / (1100*u.um)).decompose()**-x for x in range(1,5)],
                    [(3*u.mm/(1100*u.um)).decompose().value**x for x in range(1,5)], 'b--', zorder=-10)

            pl.plot(xlims, [(3*u.mm/(1100*u.um)).decompose().value**3]*2, 'k--', zorder=-10)
            pl.plot(xlims, [(3*u.mm/(1100*u.um)).decompose().value**4]*2, 'k:', zorder=-10)
            pl.plot(xlims, [(3*u.mm/(1100*u.um)).decompose().value**2]*2, 'k-', zorder=-10)
            #pl.plot([(6*u.cm / (3*u.mm)).decompose()**2]*2, [1,1e5], 'g--', zorder=-10)
            pl.plot([((350*u.um) / (1100*u.um)).decompose()**-2]*2, ylims, 'r-', zorder=-10)
            pl.plot([((350*u.um) / (1100*u.um)).decompose()**-3]*2, ylims, 'r:', zorder=-10)
            pl.plot([((350*u.um) / (1100*u.um)).decompose()**-4]*2, ylims, 'r--', zorder=-10)
            pl.plot([((350*u.um) / (1100*u.um)).decompose()**-3.5]*2, ylims, 'r-.', zorder=-10)

            pl.gca().set_xlim(*xlims)
            pl.gca().set_ylim(*ylims)
            pl.xlabel("350 $\mu$m / 1100 $\mu$m")
            pl.ylabel("1100 $\mu$m / 3 mm")
            pl.text(15, 6, "$\\beta=0$")
            pl.text(15, 66, "$\\beta=2$")
            pl.text(50, 3, "$\\beta=1.5$", rotation='vertical', color='r')

            fign = f'{catalog_figure_path}/colorcolor_350um1100um_vs_1100um3mm_{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.pdf'


            break
        break
    break
raise

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
