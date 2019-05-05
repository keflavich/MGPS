import numpy as np
from astropy.table import Table
from paths import catalog_path, catalog_figure_path

import pylab as pl


full_table = Table.read(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac')

bolocamdetected = ~((full_table['Fint1100um'] == 0))
atlasgaldetected = ~((full_table['Fint870um'] == 0))
higaldetected = full_table['HerschelDetected'] == 'True'
cmdetected = full_table['cmDetected'] == 'True'
mmdetected = bolocamdetected | atlasgaldetected | higaldetected
cm_mm_nondetection = (~cmdetected) & (~mmdetected)


pl.clf()
bins = np.logspace(np.log10(2.5e-3), 1)
pl.hist(full_table['MUSTANG_dend_flux'], bins=bins, log=True,
        label="All sources")
pl.hist(full_table['MUSTANG_dend_flux'][cm_mm_nondetection],
        bins=bins, log=True, label="cm/mm nondetections")
pl.semilogx()
pl.legend(loc='best')
pl.xlabel("MUSTANG source flux $S_{3 \mathrm{mm}}$ [Jy]")
pl.savefig(f'{catalog_figure_path}/full_catalog_histogram.pdf')
