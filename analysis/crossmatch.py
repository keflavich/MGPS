import numpy as np
from astropy.table import Table, Column
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy import wcs

from paths import catalog_path
from files import files
import constants

higal_catalog = 'http://vialactea.iaps.inaf.it/vialactea/public/HiGAL_clump_catalogue_v1.tar.gz'
catalogs_to_search = {'J/AJ/131/2525/table2': {'Fpeak':'Fpeak20cm', #mJy
                                               'Fint':'Fint20cm',
                                               'RMS':'RMS20cm'},
                      'J/A+A/619/A124/catalog': {'Sp': 'Fpeak20cm_THOR', # Jy/beam
                                                 'Sint': 'Fint20cm_THOR',
                                                },
                      'J/AJ/129/348/cat6cm': {'F5GHzp': 'Fpeak6cm_MAGPIS', }, # mJy
                      'J/ApJS/205/1/catalog': {'Sint': 'Fint6cm_CORNISH',
                                               'e_Sint': 'eFint6cm_CORNISH',
                                              },
                      'J/ApJS/91/347/table2': {'Sp5GHz': 'Fpeak6cm_Becker',
                                               #'Si5GHz': 'Fint6cm_Becker',
                                              },
                      'J/A+A/568/A41/atl-csc': {'Fint': 'Fint870um',
                                                'e_Fint': 'e_Fint870um',},
                      #'J/MNRAS/473/1059/table5': {'Fint': 'Fint870um', # Jy
                      #                            'Fpeak': 'Fpeak870um', # Jy/beam
                      #                           }
                      #'J/ApJ/799/29/table5': {
                      'J/AJ/149/64/catalog': {'S24': 'Fint24um', # mJy
                                              '__8.0_': 'Fint8um', #mag
                                              '__5.8_': 'Fint5_8um',
                                              '__4.5_': 'Fint4_5um',
                                              '__3.6_': 'Fint3_6um',
                                             },
                      'J/ApJS/208/11/table1': {'Type': 'RMSClass'}, # RMS Lumsden+ 2013
                      # int: Jy, peak: MJy/sr
                      'J/A+A/591/A149/higalblu': {'Fint': 'Fint70um', 'Fpeak': 'Fpeak70um'}, # Blue (PACS 70um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalred': {'Fint': 'Fint160um', 'Fpeak': 'Fpeak160um'}, # Red (PACS 160um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpsw': {'Fint': 'Fint250um', 'Fpeak': 'Fpeak250um'}, # PSW (SPIRE 250um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpmw': {'Fint': 'Fint350um', 'Fpeak': 'Fpeak350um'}, # PMW (SPIRE 350um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalplw': {'Fint': 'Fint500um', 'Fpeak': 'Fpeak500um'}, # PMW (SPIRE 350um) band HIGAL Herschel catalog
                      'IRSA:bolocamv21': {'flux': 'Fint1100um', 'flux_40': 'Fint1100um_40as', 'err_flux_40': 'eFint1100um_40as'}, # Jy
                      'J/A+A/588/A97/catalog': {'alpha': 'alpha_THOR', 'e_alpha': 'e_alpha_THOR',
                                                'Sp': 'Fpeak20cm_THOR2',} # Jy
                     }
#magpis_both_catalog = 'J/AJ/130/586'

if __name__ == "__main__":
    for regname,fn in files.items():
        print(f"Crossmatching region {regname}, file {fn}")
        header = fits.getheader(fn)
        ww = wcs.WCS(header)
        frame = wcs.utils.wcs_to_celestial_frame(ww)

        # loop over each catalog (we ended up with one, but early on we had
        # others to test with)
        for threshold,min_npix in ((4, 100),):# (4, 15)): #(6, 15), (8, 15), (10, 15)):
            for min_delta in (1, ): #2):

                ppcat = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.ipac', format='ascii.ipac')

                # for each table, we need to add the column from the search
                # target.  The RMS catalog gets special-cased because the first
                # attempt resulted in string lengths that were too short
                for vcatname, coldesc in catalogs_to_search.items():
                    for colname in coldesc.values():
                        if 'RMSClass' in colname:
                            ppcat.add_column(Column(name=colname, length=len(ppcat), dtype='S20'))
                        else:
                            ppcat.add_column(Column(name=colname, length=len(ppcat)))

                # for each non-rejected data point in the catalog, do the searches
                for row in ppcat:
                    if row['rejected'] == 0 or row['rejected'] == 'False':
                        crd = coordinates.SkyCoord(row['x_cen'], row['y_cen'],
                                                   frame=frame.name,
                                                   unit=(u.deg, u.deg))

                        # for each catalog we're going to search, run the query
                        for vcat,coldesc in catalogs_to_search.items():

                            # IRSA and Vizier are the two cases we're searching
                            if 'IRSA' in vcat:
                                rslt = Irsa.query_region(crd,
                                                         radius=10*u.arcsec,
                                                         catalog=vcat.split(":")[1])
                                if len(rslt) == 1:
                                    rslt = [rslt] # hack b/c Irsa returns 1 table, vizier returns a list
                            else:
                                rslt = Vizier.query_region(crd,
                                                           radius=10*u.arcsec,
                                                           catalog=vcat)

                            #this block does nothing: we already just cut to the 0'th
                            # # if there are multiple hits, we don't have any way
                            # # to decide between them, so we just skip it
                            # # (this isn't great; ideally we'd like to do something about this...)
                            # if len(rslt) > 1:
                            #     print(f"WARNING: cropping result down from {rslt} to {rslt[:1]}")
                            #     rslt = rslt[:1]
                            #     #pass
                            #     #print(rslt)

                            if len(rslt) == 1:
                                # convert from tabledict to table
                                tbl = rslt[0]
                                tblrow = tbl[0]
                                for origcolname,colname in coldesc.items():
                                    row[colname] = tblrow[origcolname]

                                    # make sure we get the units set correctly
                                    if ppcat[colname].unit is None and tbl[origcolname].unit is not None:
                                        ppcat[colname].unit = tbl[origcolname].unit
                                        print(f"Set unit for {colname} to {tbl[origcolname].unit}")
                                    if ('alpha' in origcolname) or ('Type' in origcolname):
                                        pass
                                    elif tbl[origcolname].unit is None:
                                        raise

                                    # sanity check; turns out we can't reach
                                    # this because there are no hits for some
                                    # fields (not all fields overlap w/the
                                    # survey)
                                    if 'Fint20cm' in (colname, origcolname):
                                        assert tbl[origcolname].unit is not None
                                        assert ppcat[colname].unit is not None

                herscheldetected = (ppcat['Fint70um', 'Fint160um', 'Fint250um', 'Fint350um'].as_array().view('float').reshape(len(ppcat),4) > 0).any(axis=1)
                spitzerdetected = (ppcat['Fint8um', 'Fint3_6um', 'Fint4_5um', 'Fint5_8um', 'Fint24um'].as_array().view('float').reshape(len(ppcat),5) > 0).any(axis=1)
                cmdetected = (ppcat['Fint6cm_CORNISH', 'Fpeak6cm_Becker', 'Fpeak6cm_MAGPIS','Fint20cm', 'Fint20cm_THOR'].as_array().view('float').reshape(len(ppcat),5) > 0).any(axis=1)
                ppcat.add_column(Column(name='HerschelDetected', data=herscheldetected))
                ppcat.add_column(Column(name='SpitzerDetected', data=spitzerdetected))
                ppcat.add_column(Column(name='cmDetected', data=cmdetected))

                # special case: what if there are _no_ hits?  Unfortunately this might have to be done for multiple columns across the different fields
                if ppcat['Fint20cm'].unit is None:
                    ppcat['Fint20cm'].unit = u.Jy
                if ppcat['Fint6cm_CORNISH'].unit is None:
                    ppcat['Fint6cm_CORNISH'].unit = u.Jy
                #if ppcat['Fint6cm_Becker'].unit is None:
                #    ppcat['Fint6cm_Becker'].unit = u.mJy
                if ppcat['Fpeak6cm_Becker'].unit is None:
                    ppcat['Fpeak6cm_Becker'].unit = u.mJy
                if ppcat['Fpeak6cm_MAGPIS'].unit is None:
                    ppcat['Fpeak6cm_MAGPIS'].unit = u.Jy
                for key in (70,160,250,350,500):
                    if ppcat[f'Fpeak{key}um'].unit is None:
                        ppcat[f'Fpeak{key}um'].unit = u.Jy/u.sr

                ppcat['x_cen'].unit = u.deg
                ppcat['y_cen'].unit = u.deg

                ppcat['3mm20cmindex_THOR'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fpeak20cm_THOR'])) / np.log(constants.mustang_central_frequency / (20*u.cm).to(u.GHz, u.spectral()))
                ppcat['3mm20cmindex'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fpeak20cm'])) / np.log(constants.mustang_central_frequency / (20*u.cm).to(u.GHz, u.spectral()))
                ppcat['3mm1mmindex'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fint1100um_40as'] * 1.46)) / np.log(constants.mustang_central_frequency / (271.1*u.GHz))
                ppcat['3mm6cmindex_MAGPIS'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fpeak6cm_MAGPIS'].quantity.to(u.Jy).value)) / np.log(constants.mustang_central_frequency / (5*u.GHz))
                ppcat['3mm6cmindex_CORNISH'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fint6cm_CORNISH'].quantity.to(u.Jy).value)) / np.log(constants.mustang_central_frequency / (5*u.GHz))
                ppcat['3mm6cmindex_Becker'] = np.log(ppcat['MUSTANG_10as_peak'] / (ppcat['Fpeak6cm_Becker'].quantity.to(u.Jy).value)) / np.log(constants.mustang_central_frequency / (5*u.GHz))

                # identify HCHII candidates from criteria:
                # S_3mm > S_6cm and/or S_20cm, or nondetections at long wavelengths plus an excess over extrapolation from 1mm at beta=3
                # Dust-detected (but does not need to be a point source)
                candidate_hchii = (
                                   ((ppcat['rejected'] == 0) | (ppcat['rejected'] == 'False')) &
                                   (
                                   ( # cm slope > -0.1
                                    ((ppcat['MUSTANG_10as_peak'] > 1.75*ppcat['Fpeak6cm_MAGPIS']) & (ppcat['Fpeak6cm_MAGPIS'] > 0)) |
                                    ((ppcat['MUSTANG_10as_peak'] > 1.75*ppcat['Fint6cm_CORNISH']) & (ppcat['Fint6cm_CORNISH'] > 0)) |
                                    ((ppcat['MUSTANG_10as_peak'] > 1.75*ppcat['Fpeak6cm_Becker']) & (ppcat['Fpeak6cm_Becker'] > 0)) |
                                    ((ppcat['MUSTANG_10as_peak'] > 43*ppcat['Fpeak20cm_THOR']) & (ppcat['Fpeak20cm_THOR'] > 0)) |
                                    ((ppcat['MUSTANG_10as_peak'] > 43*ppcat['Fpeak20cm']) & (ppcat['Fpeak20cm'] > 0))
                                   ) | 
                                   # 3mm excess over 1mm extrapolation, but no 6cm/20cm detection
                                   (
                                       (ppcat['Fpeak20cm'] == 0) & (ppcat['Fpeak20cm_THOR'] == 0) & (ppcat['Fpeak6cm_MAGPIS'] == 0) &
                                   # 1.46 is the 40as -> Gaussian correction factor
                                   # 3.0 is the beta=1 case for dust (shallow)
                                   # This criterion looks for an excess at 3mm over pure dust
                                    (
                                     (ppcat['MUSTANG_10as_peak'] / (ppcat['Fint1100um_40as'] * 1.46) > (constants.mustang_central_frequency / (271.1*u.GHz))**(3.0)) &
                                     # if no 1 mm detection, assume no dust - HCHII unlikely w/o dust
                                     (ppcat['Fint1100um_40as'] != 0)
                                    )
                                   ) |
                                   (
                                   # 3mm excess over 870um extrapolation, but no 6cm/20cm detection
                                       (ppcat['Fpeak20cm'] == 0) & (ppcat['Fpeak20cm_THOR'] == 0) & (ppcat['Fpeak6cm_MAGPIS'] == 0) &
                                   # 3.0 is the beta=1 case for dust (shallow)
                                   # This criterion looks for an excess at 3mm over pure dust
                                    (
                                     (ppcat['MUSTANG_10as_peak'] / (ppcat['Fint870um']) > (constants.mustang_central_frequency / (350*u.GHz))**(3.0)) &
                                     # if no 870 um detection, assume no dust - HCHII unlikely w/o dust
                                     (ppcat['Fint870um'] != 0)
                                    )
                                   )
                                   )
                                  )
                ppcat.add_column(Column(name="HCHII_candidate", data=candidate_hchii))

                outfn = f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.ipac'
                ppcat.write(outfn, format='ascii.ipac')
                print(f"Completed file {outfn}")


    #ds9 W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr10_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr10_minn15_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr4_minn20_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr4_minn20_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr6_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr6_minn15_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr8_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr8_minn15_mind2.reg
