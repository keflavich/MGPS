import os
from dendrocat.aperture import Circle, Annulus
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy import units as u
from astropy import coordinates
from astropy.table import Column, Table, join
import regions
import pylab as pl
from paths import catalog_figure_path, catalog_path, overview_figure_path
from files import files
from gaussfit_catalog import gaussfit_catalog, gaussfits_to_table

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

            catalog = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.ipac', format='ascii.ipac')

            centers = coordinates.SkyCoord(catalog['x_cen'], catalog['y_cen'], frame='galactic', unit=(u.deg, u.deg))
            reglist = [regions.PointSkyRegion(x, meta={'text':row['SourceName']}) for x,row in zip(centers, catalog) if row['rejected'] == 0]

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
            header['BMAJ'] = 9/3600.
            header['BMIN'] = 9/3600.
            header['BPA'] = 0.
            fh.writeto(fn, overwrite=True)

            diagnostics_dir = '/Volumes/external/mgps/gaussfit_diagnostics/'
            prefix = regname+"_"
            if not os.path.exists(diagnostics_dir):
                diagnostics_dir = None
                prefix = ''

            gfit_dat = gaussfit_catalog(fn, reglist, radius=30*u.arcsec,
                                        savepath=diagnostics_dir,
                                        prefix=prefix,
                                       )

            gfit_tbl = gaussfits_to_table(gfit_dat)

            # IPAC compatibility
            gfit_tbl.rename_column("chi2/n", "chi2_n")
            gfit_tbl.rename_column("Name", "SourceName")

            merge_tbl = join(catalog, gfit_tbl, join_type='left', keys='SourceName')

            merge_tbl.write(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch_gaussfits.ipac',
                            format='ascii.ipac', overwrite=True)
