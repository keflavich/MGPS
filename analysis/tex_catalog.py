import paths
import numpy as np
from astropy.table import Table,Column
from astropy import coordinates
from astropy import units as u
from latex_info import (latexdict, format_float, round_to_n, rounded,
                        rounded_arr, strip_trailing_zeros, exp_to_tex)
from astropy.io import fits
from astropy import wcs

latexdict = latexdict.copy()

orig_cont_tbl = Table.read(paths.tpath("concatenated_catalog.ipac"), format='ascii.ipac')
cont_tbl = Table.read(paths.tpath("concatenated_catalog.ipac"), format='ascii.ipac')

cont_tbl = cont_tbl[(cont_tbl['rejected'] == 0)]

description = {'_idx': 'Source ID number',
               'MUSTANG_10as_sum': 'MUSTANG aperture flux within a 10" aperture',
               'MUSTANG_15as_sum': 'MUSTANG aperture flux within a 15" aperture',
               'x_cen': 'Dendrogram moment-1 Galactic longitude',
               'y_cen': 'Dendrogram moment-1 Galactic latitude',
               'MUSTANG_dend_flux': 'Dendrogram contour-integrated intensity',
               'MUSTANG_background_median': 'Median background intensity in 15"-20" annulus',
               'amplitude': 'Gaussian fit amplitude',
               'center_x': 'Gaussian fit Galactic longitude',
               'center_y': 'Gaussian fit Galactic latitude',
               'fwhm_major': 'Gaussian fit FWHM along the major axis',
               'fwhm_minor': 'Gaussian fit FWHM along the minor axis',
               'pa': 'Gaussian fit position angle',
               'e_amplitude': 'error on Gaussian fit amplitude',
               'e_center_x': 'error on Gaussian fit Galactic longitude',
               'e_center_y': 'error on Gaussian fit Galactic latitude',
               'e_fwhm_major': 'error on Gaussian fit FWHM along the major axis',
               'e_fwhm_minor': 'error on Gaussian fit FWHM along the minor axis',
               'e_pa': 'error on Gaussian fit position angle',
              }

cont_tbl['MUSTANG_10as_sum'].unit = u.Jy
cont_tbl['MUSTANG_15as_sum'].unit = u.Jy
cont_tbl['MUSTANG_background_median'].unit = u.Jy/u.beam
cont_tbl['amplitude'].unit = u.Jy/u.beam


rename_mapping = {'_idx': 'ID',
                  'MUSTANG_dend_flux': 'Dendrogram $S_{\\nu}$',
                  'x_cen': '$\ell$',
                  'y_cen': '$b$',
                  'MUSTANG_10as_sum': '$S_{\\nu,10\'\'}$',
                  'MUSTANG_15as_sum': '$S_{\\nu,15\'\'}$',
                  'MUSTANG_background_median': '$S_{bg;15-20\'\'}$',
                  'amplitude': '$A_G$',
                  'center_x': '$\ell_G$',
                  'center_y': '$b_G$',
                  'fwhm_major': 'FWHM$_{maj,G}$',
                  'fwhm_minor': 'FWHM$_{min,G}$',
                  'pa': 'PA$_G$',
                  'e_amplitude': '$e_{A,G}$',
                  'e_center_x': '$e_{\ell,G}$',
                  'e_center_y': '$e_{b,G}$',
                  'e_fwhm_major': '$e_{\mathrm{FWHM},maj,G}',
                  'e_fwhm_minor': '$e_{\mathrm{FWHM},min,G}',
                  'e_pa': '$e_{\mathrm{PA},G}$',
                 }

# round the data
for key in rename_mapping.keys():
    ekey = f'e_{key}'
    if ekey in cont_tbl.colnames:
        unit = cont_tbl[key].unit
        cont_tbl[key] = rounded_arr(cont_tbl[key], cont_tbl[ekey], extra=0)
        cont_tbl[key].unit = unit
        cont_tbl[ekey] = rounded_arr(cont_tbl[ekey], cont_tbl[ekey], extra=1)
        cont_tbl.remove_column(ekey)

for old, new in rename_mapping.items():
    if old in cont_tbl.colnames:
        cont_tbl[old].meta['description'] = description[old]
        cont_tbl.rename_column(old, new)


# remove all not-renamed columns
for colname in ['_idx', '_index', '_name', 'area_ellipse', 'area_exact',
                'MUSTANG_dend_flux', 'major_fwhm', 'minor_fwhm',
                'position_angle', 'radius', 'x_cen', 'y_cen', 'rejected',
                'MUSTANG_detected', 'MUSTANG_snr', 'MUSTANG_10as_peak',
                'MUSTANG_10as_sum', 'MUSTANG_10as_rms', 'MUSTANG_10as_median',
                'MUSTANG_10as_npix', 'MUSTANG_15as_peak', 'MUSTANG_15as_sum',
                'MUSTANG_15as_rms', 'MUSTANG_15as_median', 'MUSTANG_15as_npix',
                'MUSTANG_background_peak', 'MUSTANG_background_sum',
                'MUSTANG_background_rms', 'MUSTANG_background_median',
                'MUSTANG_background_npix', 'SourceName', 'Fpeak20cm',
                'Fint20cm', 'RMS20cm', 'Fpeak20cm_THOR', 'Fint20cm_THOR',
                'Fpeak6cm_MAGPIS', 'Fint6cm_CORNISH', 'eFint6cm_CORNISH',
                'Fint6cm_CORNISH_2', 'CORNISH_EM', 'CORNISH_tau', 'CORNISH_ne',
                'Fpeak6cm_Becker', 'Fint870um', 'e_Fint870um', 'Fint24um',
                'Fint8um', 'Fint5_8um', 'Fint4_5um', 'Fint3_6um', 'RMSClass',
                'Fint70um', 'Fpeak70um', 'Fint160um', 'Fpeak160um',
                'Fint250um', 'Fpeak250um', 'Fint350um', 'Fpeak350um',
                'Fint500um', 'Fpeak500um', 'Fint1100um', 'Fint1100um_40as',
                'eFint1100um_40as', 'alpha_THOR', 'e_alpha_THOR',
                'Fpeak20cm_THOR2', 'HerschelDetected', 'SpitzerDetected',
                'cmDetected', '3mm20cmindex_THOR', '3mm20cmindex',
                '3mm1mmindex', '3mm6cmindex_MAGPIS', '3mm6cmindex_CORNISH',
                '3mm6cmindex_Becker', 'HCHII_candidate', 'amplitude',
                'center_x', 'center_y', 'fwhm_major', 'fwhm_minor', 'pa',
                'chi2', 'chi2_n', 'e_amplitude', 'e_center_x', 'e_center_y',
                'e_fwhm_major', 'e_fwhm_minor', 'e_pa', 'success',
                'AspectRatio', 'MorphologyClass', 'FieldID',
                'bolocam_and_mgps_detected', 'atlasgal_and_mgps_detected',
                'higal_and_mgps_detected', 'cm_and_mgps_detected',
                'mm_and_mgps_detected', 'mgps_detected_mm_and_cm_nondetected',
                'mgps_and_mm_detected_cm_nondetected',
                'mgps_and_mm_detected_nocm_compact']:
    if colname in cont_tbl.colnames:
        cont_tbl.remove_column(colname)


formats = {key: lambda x: ('{0:0.2f}'.format(np.round(x,2)))
           for key in rename_mapping.values()}

formats.update({#'Coordinates': lambda x: x.to_string('hmsdms', sep=":"),
           '$S_{\\nu,10\'\'}$': lambda x: strip_trailing_zeros(str(x)), #'{0:0.2f}'.format(round_to_n(x,2))),
           '$S_{\\nu,15\'\'}$': lambda x: strip_trailing_zeros(str(x)), #'{0:0.2f}'.format(round_to_n(x,2))),
           '$S_{bg;15-20\'\'}$': lambda x: strip_trailing_zeros(str(x)), #'{0:0.2f}'.format(round_to_n(x,2))),
           '$\ell$': lambda x: ('{0:0.3f}'.format(np.round(x,3))),
           '$b$': lambda x: ('{0:0.3f}'.format(np.round(x,3))),
           '$\ell_G$': lambda x: ('{0:0.3f}'.format(np.round(x,3))),
           '$b_G$': lambda x: ('{0:0.3f}'.format(np.round(x,3))),
           'PA$_G$': lambda x: strip_trailing_zeros(str(x)),# lambda x: ('{0:0.1f}'.format(np.round(x,1))),
           'FWHM$_{maj,G}$': lambda x: str(x),# lambda x: ('{0:0.1f}'.format(np.round(x,1))),
           'FWHM$_{min,G}$': lambda x: str(x),# lambda x: ('{0:0.1f}'.format(np.round(x,1))),
          })

# shorter-form units
#cont_tbl['$S_{\\nu,max}$'].unit = 'mJy bm$^{-1}$'
#cont_tbl['$\sigma_{bg}$'].unit = 'mJy bm$^{-1}$'


for col in formats:
    #if col == 'Coordinates':
    #    continue
    if col not in cont_tbl.colnames:
        continue
    if hasattr(cont_tbl[col], 'unit') and cont_tbl[col].unit is not None:
        cont_tbl[col] = u.Quantity(list(map(lambda x: np.round(x, 3), cont_tbl[col])),
                                   cont_tbl[col].unit)
        assert hasattr(cont_tbl[col], 'unit')
        assert cont_tbl[col].unit is not None
    else:
        cont_tbl[col] = list(map(lambda x: np.round(x, 3), cont_tbl[col]))


cont_tbl.write(paths.tpath('continuum_photometry_forpub.tsv'), format='ascii.csv', delimiter='\t', overwrite=True)

# coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'])
# cont_tbl.remove_column('RA')
# cont_tbl.remove_column('Dec')
# cont_tbl.add_column(Column(name="Coordinates", data=coords))
# cont_tbl.remove_column('SIMBAD_ID')


# caption needs to be *before* preamble.
#latexdict['caption'] = 'Continuum Source IDs and photometry'
latexdict['header_start'] = '\label{tab:photometry}'#\n\\footnotesize'
latexdict['preamble'] = '\caption{\\MUSTANG Source IDs and photometry}\n\\resizebox{\\textwidth}{!}{'
latexdict['col_align'] = 'l'*len(cont_tbl.columns)
latexdict['tabletype'] = 'table'
latexdict['tablefoot'] = ("}\par\n"
                          "The subscripts X${_G}$ are for the parameters derived"
                          " from Gaussian fits.  The values displayed are rounded "
                          "such that the error is in the last digit; error estimates "
                          "can be found in the digital version of the table."
                          "Note that position angles in the set (0, 90, 180, 270, 360) "
                          "are caused by bad fits.  These fits are kept in the catalog "
                          "because they passed other criteria and are high signal-to-noise, "
                          "but they are likely of sources in crowded regions so the "
                          "corresponding fit parameters should be treated with caution."

                         )

cont_tbl.sort('$S_{\\nu,10\'\'}$')
cont_tbl[:-20:-1].write(paths.texpath("continuum_photometry.tex"),
                        formats=formats, overwrite=True, latexdict=latexdict)

