import os

rootdir = '/Users/adam/work/mgps/'

figure_path = rootdir + 'figures/'
overview_figure_path = rootdir + 'figures/overview/'
catalog_figure_path = rootdir + 'figures/cataloging/'
catalog_path = rootdir + 'tables/'
basepath = rootdir + ''
pilotpaperpath = rootdir + 'pilotpaper'


def root(x, rootdir=rootdir):
    return os.path.join(rootdir, x)
