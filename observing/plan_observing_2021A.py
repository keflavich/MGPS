import numpy as np
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

assert len(targets) < 65

# observability tables take _forever_
# first do a fast one to make sure things work...
#obstab_quick = observability_table(constraints, observer, targets[0:2],
#                                   time_range=time_range,
#                                   time_grid_resolution=1*u.hour)
#print(obstab_quick)
#
#
#observability_table = observability_table(constraints, observer, targets[2:],
#                                          time_range=time_range,
#                                          time_grid_resolution=1*u.hour)
#print(observability_table)

## calculate observability of most-observable target
#observability = is_observable(constraints, observer, targets[-1],
#                              time_range=time_range,
#                              time_grid_resolution=1*u.hr)


# Create grid of times from ``start_time`` to ``end_time``
# with resolution ``time_resolution``
time_resolution = 7*u.day
time_grid = time_grid_from_range(time_range, time_resolution=time_resolution)

observability_grid = np.zeros((24, len(time_grid)), dtype='bool')

for targetlon in (-5, 5, 15, 25, 35, 45, 55):
    target = coordinates.SkyCoord(targetlon*u.deg, 0*u.deg, frame='galactic')

    pb = ProgressBar(observability_grid.shape[1])

    for ii, day in (enumerate(time_grid)):
        day_time_grid = time_grid_from_range(Time([day, day+TimeDelta(1*u.day)]),
                                             time_resolution=1*u.hour)
        assert np.all(np.isfinite(day_time_grid.value))
        observability_grid[:, ii] = np.all([constraint(observer, target,
                                                       times=day_time_grid)
                                            for constraint in constraints],
                                           axis=0)
        pb.update()

    # Create plot showing observability of the target:

    extent = [-0.5, -0.5+len(time_grid), -0.5, 24.5]

    fig = pl.figure(1)
    fig.clf()
    ax = pl.gca()
    ax.imshow(observability_grid, extent=extent, cmap='gray_r')
    ax.set_aspect((6 * (1*u.day) / time_resolution).decompose().value)

    #ax.set_yticks(range(0, 3))
    #ax.set_yticklabels([c.__class__.__name__ for c in constraints])
    #ax.set_yticklabels

    #ax.set_xticks(np.arange(extent[0], extent[1]), minor=False)
    step = 2
    ax.set_xticks(range(0, len(time_grid), step))
    ax.set_xticklabels([t.datetime.strftime("%Y/%m/%d") for t in time_grid[::step]])

    ax.set_yticks(np.arange(extent[2], extent[3]), minor=True)

    #ax.grid(which='minor', color='w', linestyle='-', linewidth=0.5)
    #ax.tick_params(axis='x', which='minor', bottom='off')
    pl.setp(ax.get_xticklabels(), rotation=30, ha='right')

    ax.tick_params(axis='y', which='minor', left='off')
    ax.set_xlabel('Time')
    fig.subplots_adjust(left=0.35, right=0.9, top=0.9, bottom=0.1)
    #pl.show()
    pl.savefig(f"ObservabilityOfG{targetlon}.png", bbox_inches='tight')
