import numpy as np
from astropy.table import Table
from paths import catalog_path, catalog_figure_path, pilotpaperpath
import powerlaw
import plfit

import pylab as pl


full_table = Table.read(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac')

bolocamdetected = ~((full_table['Fint1100um'] == 0))
atlasgaldetected = ~((full_table['Fint870um'] == 0))
higaldetected = full_table['HerschelDetected'] == 'True'
cmdetected = full_table['cmDetected'] == 'True'
mmdetected = bolocamdetected | atlasgaldetected | higaldetected
cm_mm_nondetection = (~cmdetected) & (~mmdetected)
compact = full_table['MorphologyClass'] == 'C'


pl.figure(1)
pl.clf()
pl.plot(full_table['MUSTANG_dend_flux'] / full_table['MUSTANG_10as_peak'], full_table['fwhm_major'], '.')
pl.plot(full_table['MUSTANG_15as_sum'] / full_table['MUSTANG_15as_peak'], full_table['fwhm_major'], '.')


pl.legend(loc='best')
pl.xlabel("'Y-factor' ($S_{int}/S_{peak}$)")
pl.ylabel("Major axis FWHM (arcsec)")
pl.savefig(f'{catalog_figure_path}/full_catalog_yfactor_vs_fwhm.pdf')

pl.figure(2)
pl.clf()
pl.plot(full_table['MUSTANG_10as_peak'], full_table['MUSTANG_15as_peak'], '.')


pl.xlabel("Peak (10as)")
pl.ylabel("Peak (15as)")
pl.savefig(f'{catalog_figure_path}/full_catalog_peak_estimates.pdf')

pl.figure(3)
pl.clf()
pl.plot([0,7],[0,7],'k-')
pl.plot(full_table['MUSTANG_10as_peak'], full_table['MUSTANG_10as_sum'], '.')
pl.plot(full_table['MUSTANG_15as_peak'], full_table['MUSTANG_15as_sum'], '.')


pl.legend(loc='best')
pl.xlabel("Peak")
pl.ylabel("Sum")
pl.savefig(f'{catalog_figure_path}/full_catalog_peak_vs_sum.pdf')
