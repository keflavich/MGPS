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

            catalog = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.ipac', format='ascii.ipac')


