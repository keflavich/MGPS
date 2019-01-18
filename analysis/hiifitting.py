import pylab as pl
import imp
sys.path.append("/Users/adam/repos/hiimodel/hiimodel/")
import HII_model
imp.reload(HII_model)
from astropy import table
from astropy import units as u
from utils import makesed

ppcat = table.Table.read('../tables/G31_dend_contour_thr4_minn20_mind1_crossmatch.ipac', format='ascii.ipac')

wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.666-0.332'][0])

# this doesn't really give a fit at all
mod = HII_model.HIIregion(nu=wavelength.to(u.GHz, u.spectral()), flux=flux,
                          fluxerr=flux*0.2, normfac=1e-10,
                          alpha=4.00,
                          normfac2=1e-53, quiet=0, dust=True)
print(mod.physprops())

pl.clf()
mod.loglogplot(numax=10*u.THz, dust=True)
