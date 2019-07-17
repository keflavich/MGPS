import numpy as np
import glob
from astropy.table import Table
from astropy import units as u
import pylab as pl

tables = {}
for fn in glob.glob('t1619*txt'):
    tbl = Table.read(fn, format='ascii.csv', delimiter='\t', data_start=0,
                     names=['Wavenumber','Transmission'])
    tables[fn] = tbl

xax = tables['t1619r19.txt']['Wavenumber']
filterfunc = np.product([np.interp(xax, tbl['Wavenumber'], tbl['Transmission']) for tbl in tables.values()],
                        axis=0)

frq = u.Quantity(xax, u.cm**-1).to(u.GHz, u.spectral())
filterfunc[frq < 75*u.GHz] = 0

pl.plot(frq, filterfunc)

# calculate effective central frequency
nu0 = 100*u.GHz

for alpha in np.arange(0, 4.5, 0.5):
    flux = (frq/nu0)**alpha
    ctrfrq = (frq * flux * filterfunc).sum() / (filterfunc * flux).sum()
    print(f"{alpha} & {ctrfrq:0.2f} & {ctrfrq.to(u.mm, u.spectral()):0.3f}")


"""
effective central frequencies
0.0: 89.72196996399734 GHz
0.5: 90.17848572026456 GHz
1.0: 90.63580283990281 GHz
1.5: 91.09266772612499 GHz
2.0: 91.54784678125006 GHz
2.5: 92.00014295969581 GHz
3.0: 92.44841115744173 GHz
3.5: 92.89157202413543 GHz
4.0: 93.32862386559712 GHz
"""
