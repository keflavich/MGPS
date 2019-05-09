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
            tables_.append(ppcat)


full_table = table.vstack(tables_)

full_table.write(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac')

print(f"Total of {len(full_table)} sources found.")

bolocamdetected = ~((full_table['Fint1100um'] == 0))
atlasgaldetected = ~((full_table['Fint870um'] == 0))
higaldetected = full_table['HerschelDetected'] == 'True'
cmdetected = full_table['cmDetected'] == 'True'
mmdetected = bolocamdetected | atlasgaldetected | higaldetected
cm_mm_nondetection = (~cmdetected) & (~mmdetected)
cm_no_mm_yes = (~cmdetected) & (mmdetected)


extended = full_table['fwhm_major'] > 14
compact = full_table['fwhm_major'] < 14
filamentary = full_table['fwhm_major'] / full_table['fwhm_minor'] > 1.5
labelcol = np.array(['E']*len(extended))
labelcol[compact] = 'C'
labelcol[filamentary] = 'F'
full_table.add_column(Column(name='MorphologyClass', data=labelcol))

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
    fh.write(r"\newcommand{\ncompact}{"+str(compact.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\nextended}{"+str(extended.sum())+r"\xspace}""\n")
    fh.write(r"\newcommand{\nfilamentary}{"+str(filamentary.sum())+r"\xspace}""\n")

print("Brightest nondetections:")
print(full_table[cm_mm_nondetection & (full_table['MUSTANG_dend_flux'] > 0.1)]['SourceName','RMSClass'])
