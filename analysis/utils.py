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
