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

# include a component for the Ruze equation...
def ruze(nu, epsilon=0.23025861*u.mm, eta_a=0.71):
    # 0.23 = 0.71 e^((-4 pi epsilon / 110 GHz in mm)^2) from Frayer+2018
    # epsilon = 0.23025861 mm = (-np.log(0.23 / 0.71) / (4*np.pi)**2 * (110*u.GHz).to(u.mm, u.spectral())**2) ** 0.5
    return (eta_a * np.exp(-(4*np.pi*epsilon / (nu.to(u.mm, u.spectral())))**2)).decompose()

for alpha in np.arange(0, 4.5, 0.5):
    flux = (frq/nu0)**alpha
    ruze_component = ruze(frq)
    ctrfrq = (frq * flux * filterfunc * ruze_component).sum() / (filterfunc * ruze_component * flux).sum()
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

with Ruze:
0.0 & 88.37 GHz & 3.392 mm
0.5 & 88.81 GHz & 3.376 mm
1.0 & 89.26 GHz & 3.359 mm
1.5 & 89.71 GHz & 3.342 mm
2.0 & 90.16 GHz & 3.325 mm
2.5 & 90.61 GHz & 3.309 mm
3.0 & 91.06 GHz & 3.292 mm
3.5 & 91.51 GHz & 3.276 mm
4.0 & 91.95 GHz & 3.260 mm
"""
