import numpy as np
from astropy import table
from astropy.table import Table, Column
from paths import catalog_path
from files import files



tables_ = []
for regname,fn in files.items():
    for threshold,min_npix in ((4, 100),):# (4, 15)): #(6, 15), (8, 15), (10, 15)):
        for min_delta in (1, ): #2):

            ppcat = Table.read(f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch_gaussfits.ipac', format='ascii.ipac')
            # not all tables have same length source names =(
            ppcat['SourceName'] = ppcat['SourceName'].astype('<U16')
            ppcat.add_column(Column(name='FieldID', data=[regname]*len(ppcat)))
            tables_.append(ppcat)


full_table = table.vstack(tables_)

full_table.write(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac', overwrite=True)

print(f"Total of {len(full_table)} sources found.")

mgps_ok = full_table['rejected'] == 0

bolocamdetected = ~((full_table['Fint1100um'] == 0)) & mgps_ok
atlasgaldetected = ~((full_table['Fint870um'] == 0)) & mgps_ok
higaldetected = (full_table['HerschelDetected'] == 'True') & mgps_ok
cmdetected = (full_table['cmDetected'] == 'True') & mgps_ok
mmdetected = (bolocamdetected | atlasgaldetected | higaldetected) & mgps_ok
cm_mm_nondetection = ((~cmdetected) & (~mmdetected)) & mgps_ok
cm_no_mm_yes = ((~cmdetected) & (mmdetected)) & mgps_ok


full_table.add_column(Column(name='bolocam_and_mgps_detected', data=bolocamdetected))
full_table.add_column(Column(name='atlasgal_and_mgps_detected', data=atlasgaldetected))
full_table.add_column(Column(name='higal_and_mgps_detected', data=higaldetected))
full_table.add_column(Column(name='cm_and_mgps_detected', data=cmdetected))
full_table.add_column(Column(name='mm_and_mgps_detected', data=mmdetected))
full_table.add_column(Column(name='mgps_detected_mm_and_cm_nondetected', data=cm_mm_nondetection))
full_table.add_column(Column(name='mgps_and_mm_detected_cm_nondetected', data=cm_no_mm_yes))


# extended = full_table['fwhm_major'] > 14
# compact = full_table['fwhm_major'] < 14
# filamentary = full_table['fwhm_major'] / full_table['fwhm_minor'] > 1.5
#labelcol = np.array(['E']*len(extended))
#labelcol[compact] = 'C'
#labelcol[filamentary] = 'F'
#full_table.add_column(Column(name='MorphologyClass', data=labelcol))
extended = full_table['MorphologyClass'] == 'E'
filamentary = full_table['MorphologyClass'] == 'F'
compact = full_table['MorphologyClass'] == 'C'

compact_nocm_yesmm = (~cmdetected) & (mmdetected) & compact

full_table.add_column(Column(name='mgps_and_mm_detected_nocm_compact', data=compact_nocm_yesmm))

full_table.write(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac', overwrite=True)

print(f"Bolocam: {bolocamdetected.sum()}")
print(f"ATLASGAL: {atlasgaldetected.sum()}")
print(f"higal: {higaldetected.sum()}")
print(f"cm: {cmdetected.sum()}")
print(f"mm: {mmdetected.sum()}")

with open('../pilotpaper/nsources.tex', 'w') as fh:
    fh.write(r"\newcommand{\nsources}{"+str(len(full_table))+r"\xspace}""\n")
    fh.write(r"\newcommand{\cmdetections}{"+str(cmdetected.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\mmdetections}{"+str(mmdetected.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\cmmmnondetections}{"+str(cm_mm_nondetection.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\mmdetectionscmnondetections}{"+str(cm_no_mm_yes.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\mmdetectionscmnondetectionscompact}{"+str(compact_nocm_yesmm.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\ncompact}{"+str(compact.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\nextended}{"+str(extended.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\nfilamentary}{"+str(filamentary.sum())+r"\xspace}""\n")

print("Brightest cm nondetections:")
print(full_table[cm_mm_nondetection & (full_table['MUSTANG_dend_flux'] > 0.1)]['SourceName','RMSClass','MorphologyClass'])

print("Compact cm nondetections: ")
print(full_table[compact_nocm_yesmm]['SourceName','RMSClass','MorphologyClass'])

print(" ".join([
      f"figures/cataloging/seds/SED_plot_{row['FieldID']}_{row['SourceName']}.png"
      for row in full_table[compact_nocm_yesmm]
]))


candidates = {
    'G30.348+0.392': 'YSO/Core', # http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=30.348+%2B0.392&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G30.898+0.162': 'YSO/Core', # http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%406090911&Name=%5bWWS2012%5d%20G030.90%2b00.16&submit=submit
    'G30.943+0.035': 'RMS Evolved star', # OH/IR star http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=30.943+%2B0.035&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G11.943-0.156': 'YSO/Core', # http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=11.943+-0.156&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G12.026-0.031': 'RMS YSO', #http://rms.leeds.ac.uk/cgi-bin/public/RMS_SEARCH_RESULTS.cgi?text_field_1=G012.0260-00.0317
    'G12.112-0.126':'Maser/mm core', #http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=G12.112+-0.126&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G0.711-0.039':'HII region Sgr B2 HII T', #http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=0.711+-0.039&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G0.825-0.189':'Radio source', # why wasn't this cross-matched? http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402414360&Name=GPSR5%200.824-0.190&submit=submit
    'G359.950+0.360':'imaging artifact', # NOT real.
    'G34.096+0.019':'mm core / bubble', # http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=G34.096+%2B0.019&CooFrame=Gal&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=1&Radius.unit=arcmin&submit=submit+query&CoordList=
    'G34.412+0.236':'YSO/core', # http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%402657932&Name=NAME%20G34.4MM&submit=submit
}
