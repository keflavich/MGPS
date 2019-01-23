import numpy as np
import sys
import pylab as pl
import imp
sys.path.append("/Users/adam/repos/hiimodel/hiimodel/")
import HII_model
imp.reload(HII_model)
from astropy import table
from astropy import units as u
from utils import makesed

import paths
from catalog_flux_limits import flux_limits


def modelfit(wavelength, flux, emguess=8e6, betanormfac=1e-5, hiinormfac=8e-11,
             dusttem=45, **kwargs):

    mod2 = HII_model.HIIregion(nu=wavelength.to(u.GHz, u.spectral()), flux=flux,
                               fluxerr=flux*0.2,
                               normfac=hiinormfac,
                               EMguess=emguess,
                               beta=1.75,
                               normfac2=betanormfac, quiet=0, dustT=dusttem*u.K)
    print(mod2.params)
    print(mod2.physprops())

    pl.figure(2)
    pl.clf()
    mod2.loglogplot(numax=10*u.THz, dustT=True, **kwargs)
    ylim = pl.gca().get_ylim()
    #mod2.loglogplot(numax=10*u.THz, dustT=True, params=[8e6, 8e-11, 1.75, 1e-5, 45])
    pl.plot(np.logspace(0, 5), HII_model.inufit(em=mod2.params[0], normfac=mod2.params[1])(np.logspace(0, 5)),
            color='r', linestyle=':', zorder=-5)
    pl.plot(np.logspace(0, 5),
            (HII_model.dust.blackbody.modified_blackbody(nu=np.logspace(0, 5)*u.GHz,
                                                         temperature=mod2.dustT*u.K,
                                                         column=u.Quantity(1e18, u.cm**-2),
                                                         beta=mod2.beta) *
             u.sr).to(u.mJy).value * mod2.normfac2,
            color='orange', linestyle=':', zorder=-5)

    pl.plot([x.to(u.GHz, u.spectral()).value for x,y in zip(wavelength, flux) if not np.isfinite(y)],
            [flux_limits[x].value for x,y in zip(wavelength, flux) if not np.isfinite(y)],
            linestyle='none',
            marker='v',
            color='orange')
    pl.ylim(max([ylim[0], 1]), ylim[1])


    mod = HII_model.HIIregion(nu=wavelength.to(u.GHz, u.spectral()), flux=flux,
                              fluxerr=flux*0.2, normfac=1.8e-10,
                              EMguess=emguess,
                              alpha=2.75,
                              normfac2=3e-5, quiet=0, dust=True)
    print(mod.params)
    print(mod.physprops())

    pl.figure(1)
    pl.clf()
    mod.loglogplot(numax=10*u.THz, dust=True, **kwargs)
    #mod.loglogplot(numax=10*u.THz, dust=True, params=[3e6, 1.8e-10, 2.75, 3e-5])
    pl.plot(np.logspace(0, 5), HII_model.inufit(mod.params[0], mod.params[1])(np.logspace(0, 5)),
            color='r', linestyle=':', zorder=-5)
    ylim = pl.gca().get_ylim()

    pl.plot(np.logspace(0, 5), (np.logspace(0, 5)/1.5)**mod.alpha * mod.params[3],
            color='g', linestyle=':', zorder=-5)

    pl.plot(np.logspace(0, 5), (np.logspace(0, 5)/1.5)**mod.alpha
            * mod.params[3] + HII_model.inufit(mod.params[0],
                                               mod.params[1])(np.logspace(0, 5)),
            color='orange', linestyle='--', zorder=1)

    pl.ylim(*ylim)

    return mod, mod2

if __name__ == "__main__":

    ppcat = table.Table.read('../tables/G31_dend_contour_thr4_minn20_mind1_crossmatch.ipac', format='ascii.ipac')

    if False:
        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.666-0.332'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.666-0.332.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G29.956-0.017'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G29.956-0.017.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.009-0.273'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.009-0.273.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.022+0.157'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.022+0.157.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.212+0.429'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.212+0.429.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.252+0.054'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.252+0.054.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.534+0.021'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.534+0.021.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.589-0.042'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.589-0.042.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.720-0.082'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.720-0.082.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.754-0.049'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.754-0.049.png", bbox_inches='tight')

        wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G30.866+0.114'][0])
        alphamodel, betamodel = modelfit(wavelength, flux)
        pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G30.866+0.114.png", bbox_inches='tight')

    wavelength, flux = makesed(ppcat[ppcat['SourceName'] == 'G31.213-0.181'][0])
    alphamodel, betamodel = modelfit(wavelength, flux, hiinormfac=1e-12, betanormfac=1e-7, dusttem=55, emguess=8e7,)
    pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G31.213-0.181.png", bbox_inches='tight')
    alphamodel, betamodel = modelfit(wavelength, flux, hiinormfac=1e-12,
                                     betanormfac=1e-7, dusttem=55, emguess=8e7,
                                     do_annotations=False)
    pl.figure(2).savefig(f"{paths.catalog_figure_path}/seds/SED_fit_plot_G31.213-0.181_notext.png", bbox_inches='tight')
