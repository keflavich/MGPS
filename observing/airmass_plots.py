import astroplan
from astroplan.plots import plot_altitude, plot_schedule_airmass
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling; from astropy.io import fits, ascii; from astropy.table import Table
from astroplan import Observer
import astropy.time


targets = coordinates.SkyCoord([0.67, 10.47, 29], [0.0, 0.0, 0.0], frame='galactic', unit=(u.deg, u.deg))
observer = Observer.at_site('VLA', timezone='US/Mountain')
time = astropy.time.Time('2019-01-31', location=observer.location)

import pylab as pl
pl.clf()
for target in targets:
    plot_altitude(target, observer, time, brightness_shading=True)
pl.ylim(0,55)
pl.ylabel("Altitude")
pl.title(f"{time.strftime('%D')}")
