from astropy.table import Table, Column
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa
from astropy import coordinates
from astropy import units as u

from paths import catalog_path
from files import files

higal_catalog = 'http://vialactea.iaps.inaf.it/vialactea/public/HiGAL_clump_catalogue_v1.tar.gz'
catalogs_to_search = {'J/AJ/131/2525/table2': {'Fpeak':'Fpeak20cm', #mJy
                                               'Fint':'Fint20cm',
                                               'RMS':'RMS20cm'},
                      'J/AJ/129/348/cat6cm': {'F5GHzp': 'Fpeak6cm_MAGPIS', }, # mJy
                      'J/ApJS/205/1/catalog': {'Sint': 'Fint6cm_CORNISH',
                                               'e_Sint': 'eFint6cm_CORNISH',
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
                      # int: Jy, peak: MJy/sr
                      'J/A+A/591/A149/higalblu': {'Fint': 'Fint70um', 'Fpeak': 'Fpeak70um'}, # Blue (PACS 70um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalred': {'Fint': 'Fint160um', 'Fpeak': 'Fpeak160um'}, # Red (PACS 160um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpsw': {'Fint': 'Fint250um', 'Fpeak': 'Fpeak250um'}, # PSW (SPIRE 250um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpmw': {'Fint': 'Fint350um', 'Fpeak': 'Fpeak350um'}, # PMW (SPIRE 350um) band HIGAL Herschel catalog
                      'IRSA:bolocamv21': {'flux': 'Fint1100um', 'flux_40': 'Fint1100um_40as', 'err_flux_40': 'eFint1100um_40as'}, # Jy
                     }
#magpis_both_catalog = 'J/AJ/130/586'

if __name__ == "__main__":
    for regname,fn in files.items():
        for threshold,min_npix in ((4, 20), (6, 15), (8, 15), (10, 15)):
            for min_delta in (1, 2):
                ppcat = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.ipac', format='ascii.ipac')

                for vcatname, coldesc in catalogs_to_search.items():
                    for colname in coldesc.values():
                        ppcat.add_column(Column(name=colname, length=len(ppcat)))

                for row in ppcat:
                    if row['rejected'] == 0:
                        crd = coordinates.SkyCoord(row['x_cen'], row['y_cen'],
                                                   frame='icrs', unit=(u.deg,
                                                                       u.deg))
                        for vcat,coldesc in catalogs_to_search.items():
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
                            if len(rslt) >= 1:
                                print(rslt)

                            if len(rslt) == 1:
                                # convert from tabledict to table
                                tbl = rslt[0]
                                tblrow = tbl[0]
                                for origcolname,colname in coldesc.items():
                                    row[colname] = tblrow[origcolname]
                                    if ppcat[colname].unit is None and tbl[origcolname].unit is not None:
                                        ppcat[colname].unit = tbl[origcolname].unit

                herscheldetected = (ppcat['Fint70um', 'Fint160um', 'Fint250um', 'Fint350um'].as_array().view('float').reshape(len(ppcat),4) > 0).any(axis=1)
                spitzerdetected = (ppcat['Fint8um', 'Fint3_6um', 'Fint4_5um', 'Fint5_8um', 'Fint24um'].as_array().view('float').reshape(len(ppcat),5) > 0).any(axis=1)
                cmdetected = (ppcat['Fint6cm_CORNISH', 'Fpeak6cm_MAGPIS','Fint20cm'].as_array().view('float').reshape(len(ppcat),3) > 0).any(axis=1)
                ppcat.add_column(Column(name='HerschelDetected', data=herscheldetected))
                ppcat.add_column(Column(name='SpitzerDetected', data=spitzerdetected))
                ppcat.add_column(Column(name='cmDetected', data=cmdetected))

                ppcat.write(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch.ipac', format='ascii.ipac')

    #ds9 W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr10_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr10_minn15_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr4_minn20_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr4_minn20_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr6_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr6_minn15_mind2.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr8_minn15_mind1.reg W43/GAL_031_precon_2_arcsec_pass_9.fits -region load tables/G31_dend_contour_thr8_minn15_mind2.reg
