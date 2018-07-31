import os
import glob
import numpy as np
from astropy import table
from astropy import units as u

# parse ZT2018 models

datapath = '/Volumes/external/zhang_sedmodels/yczhang04-sedfit-959e99c/Model_SEDs/model_info/'

fitstable = '/Volumes/external/zhang_sedmodels/model_info_table.fits'

if os.path.exists(fitstable):
    tbl = table.Table.read(fitstable)
else:
    alldata = []

    for fn in glob.glob('{0}/*.dat'.format(datapath)):
        x,y,z = map(int, os.path.split(fn)[-1][:8].split("_"))
        with open(fn, 'r') as fh:
            lines = fh.readlines()
            headings = lines[0].split(',')
            data = np.array(lines[1].split(), dtype='float')

        alldata.append(np.array([x,y,z] + data.tolist()))

    unit_mappings = {'msun': 'Msun', 'rsun': 'Rsun', 'lsun': 'Lsun',
                     'msun/yr': 'Msun/yr'}
    def remap_units(x):
        if x in unit_mappings:
            return unit_mappings[x]
        else:
            return x


    units=[None, None, None] + [u.Unit(remap_units(x.split()[1].strip("()"))) for x in headings]

    datacols = [table.Column(data=xx, unit=unit)
                for xx, unit in zip(np.array(alldata).T, units)]


    tbl = table.Table(data=datacols,
                      names=['xx', 'yy', 'zz'] + [x.split()[0] for x in headings],
                     )

    tbl.write(fitstable)

# determine relative lifetime of "hchii" phase, defined by:
# (a) M_star > 15 Msun
# (b) M_star < M_env
# 
# DEFINITION TWO:
# mdotd < 1e-4
hchii_times = []
for xx in range(1,16):
    for yy in range(1, 5):
        selection = (tbl['xx'] == xx) & (tbl['yy'] == yy)

        mstargrid = np.linspace(tbl['mstar'][selection].min(), tbl['mstar'][selection].max(), 100)
        menvgrid = np.linspace(tbl['massenv'][selection].min(), tbl['massenv'][selection].max(), 100)
        timegrid = np.linspace(tbl['tnow'][selection].min(), tbl['tnow'][selection].max(), 100)
        mdotgrid = np.linspace(tbl['mdotd'][selection].min(), tbl['mdotd'][selection].max(), 100)

        hchii_sel = (mstargrid > 12) & (menvgrid > mstargrid)
        hchii_sel = (mdotgrid > 1e-4) & (mstargrid > 12)
        if any(hchii_sel):
            time_in_stage = u.Quantity(timegrid[hchii_sel].max() - timegrid[hchii_sel].min(), u.yr)
            print("Time in stage {0},{1}: {2:0.1f}   ({3} points)".format(xx, yy, time_in_stage.to(u.kyr), hchii_sel.sum()))
            hchii_times += [time_in_stage] * (selection.sum())
        else:
            print("Maximum mass is {0}".format(tbl['mstar'][selection].max()))
            hchii_times += [np.nan*u.yr] * (selection.sum())

tbl.add_column(table.Column(data=u.Quantity(hchii_times), name='HCHII_time'))

# in the end, this wasn't very useful: the models are too coarsely sampled and
# predict a very wide range of HCHII times given the arbitrary criteria I
# selected.
