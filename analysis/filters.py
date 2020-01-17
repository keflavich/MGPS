import numpy as np
import glob
from astropy.table import Table
from astropy import units as u
import pandas
import pylab as pl

tables = {}
for fn in glob.glob('t1619*txt'):
    tbl = Table.read(fn, format='ascii.csv', delimiter='\t', data_start=0,
                     names=['Wavenumber','Transmission'])
    tables[fn] = tbl

xax = tables['t1619r19.txt']['Wavenumber']
filterfunc = np.product([np.interp(xax, tbl['Wavenumber'], tbl['Transmission']) for tbl in tables.values()],
                        axis=0)

# get the filters a second way, just to be sure
xf = pandas.ExcelFile('Filters_Mustang1.5_May2013.xls')
xax2 = pandas.read_excel('Filters_Mustang1.5_May2013.xls', xf.sheet_names[0])['Unnamed: 0'].data
tables2 = [pandas.read_excel('Filters_Mustang1.5_May2013.xls', sheetname)
           for sheetname in xf.sheet_names]

filterfunc2 = np.product([np.interp(xax2, tbl['Unnamed: 0'], tbl[tbl.columns[1]])
                          for tbl in tables2], axis=0)

xf2 = pandas.ExcelFile('Filter_data.xlsm')
xax3 = np.array(pandas.read_excel('Filter_data.xlsm', xf2.sheet_names[0], header=1)['#GHz'])
tables3 = [pandas.read_excel('Filter_data.xlsm', sheetname, header=1)
           for sheetname in xf2.sheet_names if sheetname[0] == 'K']

filterfunc3 = np.product([np.interp(xax3, tbl['#GHz'], tbl['transmission'])
                          for tbl in tables3], axis=0)


frq = u.Quantity(xax, u.cm**-1).to(u.GHz, u.spectral())
filterfunc[frq < 75*u.GHz] = 0
pl.figure(1).clf()
pl.plot(frq, filterfunc, label='Spreadsheet Numbers (transcribed)')

frq2 = u.Quantity(xax2, u.cm**-1).to(u.GHz, u.spectral())
filterfunc2[frq2 < 75*u.GHz] = 0
pl.plot(frq2, filterfunc2, label='Spreadsheet Numbers (xls, Simon)')

frq3 = u.Quantity(xax3, u.GHz)
filterfunc3[frq3 < 75*u.GHz] = 0
pl.plot(frq3, filterfunc3, label='Spreadsheet Numbers (xls, Charles)')


# Brian Mason's filterfunc
frq, filterfunc = np.loadtxt('m2bp.txt').T
frq = u.Quantity(frq, u.GHz)

pl.plot(frq, filterfunc, label='m2bp')
pl.plot(frq, (filterfunc)*4, label='m2bp $\\times 4$')

pl.legend(loc='upper right')
pl.xlabel("Frequency [GHz]")
pl.ylabel("Transmission Fraction")
pl.xlim(65,160)

# calculate effective central frequency
nu0 = 100*u.GHz

# include a component for the Ruze equation...
def ruze(nu, epsilon=0.23*u.mm, eta_a=0.71):
    # 0.23 = 0.71 e^((-4 pi epsilon / 110 GHz in mm)^2) from Frayer+2018
    # epsilon = 0.23025861 mm = (-np.log(0.23 / 0.71) / (4*np.pi)**2 * (110*u.GHz).to(u.mm, u.spectral())**2) ** 0.5
    return (eta_a * np.exp(-(4*np.pi*epsilon / (nu.to(u.mm, u.spectral())))**2)).decompose()

ruze_component = ruze(frq)

for ff, name in ((filterfunc, 'm2bp'), (filterfunc2, 'xls Simon'), (filterfunc3, 'xls Charles'), ):
    print()
    print(name)
    for alpha in np.arange(0, 4.5, 0.5):
        flux = (frq/nu0)**alpha
        ctrfrq = (frq * flux * filterfunc * ruze_component).sum() / (filterfunc * ruze_component * flux).sum()
        print(f"{alpha} & {ctrfrq:0.2f} & {ctrfrq.to(u.mm, u.spectral()):0.3f}\\\\")


pl.figure(2).clf()
ruze_component = ruze(frq)
pl.plot(frq, (filterfunc * ruze_component), label='m2bp')

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

Brian Mason's + Ruze:
0.0 & 87.85 GHz & 3.413 mm\\
0.5 & 88.23 GHz & 3.398 mm\\
1.0 & 88.62 GHz & 3.383 mm\\
1.5 & 89.02 GHz & 3.368 mm\\
2.0 & 89.41 GHz & 3.353 mm\\
2.5 & 89.80 GHz & 3.338 mm\\
3.0 & 90.19 GHz & 3.324 mm\\
3.5 & 90.58 GHz & 3.310 mm\\
4.0 & 90.96 GHz & 3.296 mm\\
"""
