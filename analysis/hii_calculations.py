"""
A few calculations to determine the appropriate cutoffs for HCHII candidates
"""
import numpy as np
import sys
import imp
sys.path.append("/Users/adam/repos/hiimodel/hiimodel/")
import HII_model
imp.reload(HII_model)
from astropy import units as u
from astropy import table
from astropy.table import Table, Column
from paths import catalog_path
from files import files

emission_measure = 7e7*u.cm**-6*u.pc

tau5 = HII_model.tnu(8500*u.K, 5*u.GHz, emission_measure)

print(f"Optical depth at 5 GHz with EM={emission_measure} = {tau5}")

fL,fC,fW = HII_model.inufit(em=emission_measure)([1, 5, 100]*u.GHz)
print(f"Ratio of W-band 100 GHz to L-band 1 GHz: {fW/fL} (alpha={np.log(fW/fL)/np.log(100/1)})")
print(f"Ratio of W-band 100 GHz to C-band 5 GHz: {fW/fC} (alpha={np.log(fW/fC)/np.log(100/5)})")






nhiicand = 0
ncompacthiicand = 0
candidate_names = {}
candidate_table_ = []

for regname,fn in files.items():
    for threshold,min_npix in ((4, 100),):# (4, 15)): #(6, 15), (8, 15), (10, 15)):
        for min_delta in (1, ): #2):

            tablename = f'{catalog_path}/{regname}_dend_contour_thr{threshold}_minn{min_npix}_mind{min_delta}_crossmatch_gaussfits.ipac'
            print(tablename)
            ppcat = Table.read(tablename, format='ascii.ipac')
            ppcat['SourceName'] = ppcat['SourceName'].astype('<U16')
            ppcat.add_column(Column(name='FieldID', data=[regname]*len(ppcat)))

            not_rejected = ppcat['rejected'] == 0
            compact = (ppcat['MorphologyClass'] == 'C') & (not_rejected)
            candidate_mask = (ppcat['HCHII_candidate'] == 'True') & (not_rejected)
            this_ncand = candidate_mask.sum()
            this_ncand_compact = (candidate_mask & compact).sum()
            print(f"{regname} has {this_ncand} HCHII region candidates, of which {this_ncand_compact} are compact")
            candidate_names[regname] = (ppcat['SourceName'][candidate_mask])
            if this_ncand > 0:
                candidate_table_.append(ppcat[candidate_mask])

            nhiicand += this_ncand
            ncompacthiicand += this_ncand_compact

candidate_table = table.vstack(candidate_table_)

print(f"Total of {nhiicand} HCHII region candidates")

with open('../pilotpaper/nhiicandidates.tex', 'w') as fh:
    # latex and f-strings and \n's don't mix!
    fh.write(r"\newcommand{\nhiicand}{"+str(nhiicand)+r"\xspace}")
    fh.write(r"\newcommand{\ncompacthiicand}{"+str(ncompacthiicand)+r"\xspace}")


print()
# manual IDs
cand_qual = {
#    'G30.764-0.033': 'Extended HII region',
    'G30.944+0.035': 'OH/IR star', # index is 2.99....
#    'G30.943+0.035': 'OH/IR star', # index is 2.99....
    #'G43.165-0.030': 'W49A South',
    #'G43.148+0.014': 'W49',
    #'G43.149+0.012': 'W49',
    'G43.149+0.012': 'W49',
    'G43.166-0.030': 'W49',
    'G43.167+0.010': 'W49',
#    'G49.535-0.392': 'diffuse: W51g',
#    'G49.490-0.367': 'W51 IRS 2',
#    'G49.422-0.307': 'diffuse',
#    'G48.914-0.284': 'diffuse',
#    'G0.529-0.084': 'diffuse',
#    'G0.169+0.150': 'diffuse',
    'G34.257+0.153': 'HII', # is a hypercompact HII region
#    'G49.394-0.483': 'W51e2/e8',
#    'G0.117+0.085': 'diffuse',
}

print("The following sources have not been classified by hand:")
for row in candidate_table:
    if row['SourceName'] not in cand_qual:
        print(f"{row['SourceName']}")
print()
print("If there are no sources above, that means you're done.")


print()
print(" ".join([
      f"figures/cataloging/seds/SED_plot_{row['FieldID']}_{row['SourceName'].strip()}.png"
      for row in candidate_table
     ]))
