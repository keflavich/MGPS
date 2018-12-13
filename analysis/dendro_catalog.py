import astrodendro
import dendrocat
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy import units as u
from astropy import coordinates
import regions
import pylab as pl
from paths import catalog_figure_path, catalog_path

files = {
    'G31':'../GAL_031/GAL_031_precon_2_arcsec_pass_9.fits',
    'SgrB2':'../SgrB2/SgrB2_precon_2_arcsec_pass_9.fits',
    'W33':'../W33/W33_precon_2_arcsec_cached_pass_19.fits',
    'W49':'../W49/W49_precon_2_arcsec_pass_9.fits',
    'W51':'../W51/W51_precon_2_arcsec_pass_10.fits',
}

reglist = regions.io.read_ds9('cutout_regions.reg')
cutout_regions = {reg.meta['label']: reg for reg in reglist}

def ppcat_to_regions(cat):
    center = coordinates.SkyCoord(cat['x_cen'], cat['y_cen'], unit=(u.deg, u.deg), frame='icrs')
    rad = cat['radius'].quantity
    regs = [regions.CircleSkyRegion(cen, rr) for cen, rr in zip(center, rad)]
    return regs

for regname,fn in files.items():
    # presently hard-coded to G31
    fh = fits.open(fn)
    data = fh[0].data
    header = fh[0].header
    ww = wcs.WCS(header)
    reg = cutout_regions[regname].to_pixel(ww) 
    mask = reg.to_mask()
    data_sm = convolve_fft(data, Gaussian2DKernel(15))
    data_filtered = data-data_sm
    err = mad_std(data_filtered)
    cutout = mask.multiply(data_filtered)
    ww_cutout = ww[mask.bbox.ixmin:mask.bbox.ixmax, mask.bbox.iymin:mask.bbox.iymax]

    # add beam parameters to header
    cutout_header = ww_cutout.to_header()
    cutout_header['BMAJ'] = 8/3600.
    cutout_header['BMIN'] = 8/3600.
    cutout_header['BUNIT'] = 'Jy/beam'
    cutout_header['TELESCOP'] = 'GBT'
    cutout_header['FREQ'] = '9.0e10'

    for threshold,min_npix in ((4, 20), (6, 15), (8, 15), (10, 15)):
        for min_delta in (1, 2):
            radiosource = dendrocat.RadioSource([fits.PrimaryHDU(data=cutout,
                                                                 header=cutout_header)])
            radiosource.nu = 90*u.GHz
            radiosource.freq_id = 'MUSTANG'
            radiosource.set_metadata()
            radiosource.to_dendrogram(min_value=err*threshold,
                                      min_delta=err*min_delta,
                                      min_npix=min_npix,
                                     )
            radiosource.plot_grid(skip_rejects=False,
                                  outfile=f'{catalog_figure_path}/{regname}_dendrocat_thr{threshold}_minn{min_npix}_mind{min_delta}_prerejection.png',
                                  figurekwargs={'num': 1},
                                 )
            pl.figure(1).clf()
            radiosource.autoreject(threshold=5.0)
            print("Rejected {0} sources".format(radiosource.rejected.sum()))
            radiosource.plot_grid(skip_rejects=False,
                                  outfile=f'{catalog_figure_path}/{regname}_dendrocat_thr{threshold}_minn{min_npix}_mind{min_delta}_postrejection.png',
                                  figurekwargs={'num': 1},
                                 )
            pl.figure(1).clf()
            #dend = astrodendro.Dendrogram.compute(cutout,
            #                                      min_value=err*threshold,
            #                                      min_delta=err*min_delta,
            #                                      verbose=True, min_npix=min_npix,
            #                                      wcs=ww_cutout)
            dend = radiosource.dendrogram

            pl.figure(2).clf()
            ax = pl.gca()
            ax.cla()
            pl.imshow(cutout, cmap='gray_r', interpolation='none', origin='lower',
                      vmax=0.01, vmin=-0.001)
            pltr = dend.plotter()
            for struct in dend.leaves:
                pltr.plot_contour(ax, structure=struct, colors=['r'],
                                  linewidths=[0.9], zorder=5)
                if struct.parent:
                    while struct.parent:
                        struct = struct.parent
                    pltr.plot_contour(ax, structure=struct, colors=[(0,1,0,1)],
                                      linewidths=[0.5])

            cntr = pl.gca().collections

            pl.setp([x for x in cntr if x.get_color()[0,0] == 1], linewidth=0.25)
            pl.setp([x for x in cntr if x.get_color()[0,1] == 1], linewidth=0.25)
            pl.savefig(f'{catalog_figure_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.pdf')

            # zoom only for G31
            if regname == 'G31':
                pl.axis((1125.4006254228616, 1670.3650637799306,
                         1291.6829155596627, 1871.8063499397681))
                pl.setp([x for x in cntr if x.get_color()[0,0] == 1], linewidth=0.75) # Red
                pl.setp([x for x in cntr if x.get_color()[0,1] == 1], linewidth=0.5) # Green
                pl.savefig(f'{catalog_figure_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_zoom.pdf')

            metadata = {'data_unit': u.Jy / u.beam,
                        'beam_major': 8*u.arcsec,
                        'beam_minor': 8*u.arcsec,
                        'wcs': ww_cutout,
                        'spatial_scale': wcs.utils.proj_plane_pixel_scales(ww_cutout).mean()*u.deg,
                       }
            #ppcat = astrodendro.pp_catalog(dend, metadata)
            ppcat = radiosource.to_catalog()
            ppcat.write(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.ipac', format='ascii.ipac')

            regs = ppcat_to_regions(ppcat)
            regions.write_ds9(regions=regs, filename=f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.reg')

    # only do G31 for now, since it's hard-coded
    break
