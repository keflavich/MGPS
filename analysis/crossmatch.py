from astropy.table import Table, Column
from astroquery.vizier import Vizier
from astropy import coordinates
from astropy import units as u

from paths import catalog_path
from files import files

higal_catalog = 'http://vialactea.iaps.inaf.it/vialactea/public/HiGAL_clump_catalogue_v1.tar.gz'
catalogs_to_search = {'J/AJ/131/2525/table2': {'Fpeak':'Fpeak20cm',
                                               'Fint':'Fint20cm',
                                               'RMS':'RMS20cm'},
                      'J/AJ/129/348/cat6cm': {'F5GHzp': 'Fpeak6cm_MAGPIS', },
                      'J/ApJS/205/1/catalog': {'Sint': 'Fint6cm_CORNISH',
                                               'e_Sint': 'eFint6cm_CORNISH',
                                              },
                      'J/A+A/568/A41/atl-csc': {'Fint': 'Fint870um',
                                                'e_Fint': 'e_Fint870um',},
                      #'J/ApJ/799/29/table5': {
                      'J/AJ/149/64/catalog': {'S24': 'Fint24um',
                                              '__8.0_': 'Fint8um',
                                              '__5.8_': 'Fint5.8um',
                                              '__4.5_': 'Fint4.5um',
                                              '__3.6_': 'Fint3.6um',
                                             },
                      'J/A+A/591/A149/higalblu': {'Fint': 'Fint70um', 'Fpeak': 'Fpeak70um'}, # Blue (PACS 70um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalred': {'Fint': 'Fint160um', 'Fpeak': 'Fpeak160um'}, # Red (PACS 160um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpsw': {'Fint': 'Fint250um', 'Fpeak': 'Fpeak250um'}, # PSW (SPIRE 250um) band HIGAL Herschel catalog
                      'J/A+A/591/A149/higalpmw': {'Fint': 'Fint350um', 'Fpeak': 'Fpeak350um'}, # PMW (SPIRE 350um) band HIGAL Herschel catalog
                     }
#magpis_both_catalog = 'J/AJ/130/586'

for regname,fn in files.items():
    for threshold,min_npix in ((4, 20), (6, 15), (8, 15), (10, 15)):
        for min_delta in (1, 2):
            ppcat = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}.ipac', format='ascii.ipac')

            for vcatname, coldesc in catalogs_to_search.items():
                for colname in coldesc.values():
                    ppcat.add_column(Column(name=colname, length=len(ppcat)))

            for row in ppcat:
                if ~row['rejected']:
                    for vcat,coldesc in catalogs_to_search.items():
                        rslt = Vizier.query_region(coordinates.SkyCoord(row['x_cen'],
                                                                        row['y_cen'],
                                                                        frame='icrs',
                                                                        unit=(u.deg, u.deg)),
                                                   radius=10*u.arcsec,
                                                   catalog=vcat)
                        print(rslt)

                        if len(rslt) == 1:
                            # convert from tabledict to table
                            tbl = rslt[0]
                            if len(tbl) > 1:
                                tbl = tbl[0]
                            for origcolname,colname in coldesc.items():
                                row[colname] = tbl[origcolname]
            break
