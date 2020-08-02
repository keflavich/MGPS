import numpy as np
import datetime
from astroplan.plots import plot_altitude, plot_schedule_airmass
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling
from astropy.table import Table
from astroplan import Observer
import astropy.time
from astropy.time import Time, TimeDelta
from astropy.utils.console import ProgressBar

from astroplan import time_grid_from_range
from astroplan import (AltitudeConstraint, AtNightConstraint, Constraint)
from astroplan import (is_observable, is_always_observable, months_observable,
                       observability_table)

import pylab as pl
pl.ion()

observer = Observer.at_site('GBT', timezone='US/Eastern')

time1 = astropy.time.Time('2021-02-01', location=observer.location)
time2 = astropy.time.Time('2021-07-31', location=observer.location)
time_range = astropy.time.Time([time1,time2])



GBTNight = AtNightConstraint()
GBTNight.max_solar_altitude = 5*u.deg

class GBT3hoursAfterSunset(Constraint):
    """
    Constraint requiring 3hrs setup for MUSTANG
    """
    def compute_constraint(self, times, observer, targets):

        # we want the time since the sun went below 5 deg,
        # which is enough to start cooling the dish
        sunset = observer.sun_set_time(times, horizon=5*u.deg)

        mask = times > sunset + astropy.time.TimeDelta(3*u.hour)

        return mask

constraints = [AltitudeConstraint(20*u.deg, 80*u.deg),
               GBTNight,
               GBT3hoursAfterSunset()
              ]


targets = [coordinates.SkyCoord(glon, 0, frame='galactic', unit='deg') for
           glon in range(-5, 56)]

for target in targets:
    risetime = observer.target_rise_time(time1, target, which='next', horizon=15*u.deg)
    settime = observer.target_set_time(time1, target, which='next', horizon=15*u.deg)
    lstrise = observer.local_sidereal_time(risetime)
    lstset = observer.local_sidereal_time(settime)
    print(f"l={target.l:12s}   lst rise: {lstrise.to_string(sep=':'):15s}  lst set: {lstset.to_string(sep=':'):15s}")

