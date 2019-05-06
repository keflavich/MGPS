import json
from astropy import table
from paths import catalog_path, pilotpaperpath


noisefilepath = f"{catalog_path}/noise_levels.json"
with open(noisefilepath, 'r') as fh:
    noise_levels = json.load(fh)


obstable = {
    "W51":    {'altname':'G49', 'time':  1.02, 'sessions': "01"},
    "W49":    {'altname':'G43', 'time':  1.02, 'sessions': "01, 02"},
    "SgrB2":  {'altname':'G01', 'time':  1.35, 'sessions': "02, 03, 04, 05"},
    "GAL031": {'altname':'G31', 'time':  1.49, 'sessions': "02, 03"},
    "W33":    {'altname':'G12', 'time':  1.02, 'sessions': "03"},
    "GAL029": {'altname':'G29', 'time':  1.30, 'sessions': "04,05"},
    "GAL034": {'altname':'G34', 'time':   0.5, 'sessions': "05"},
}

for key in obstable:
    obstable[key]['noiselevel'] = noise_levels[obstable[key]['altname']]['noise']


keys_sorted = sorted(obstable, key=lambda x: obstable[x]['altname'])

with open(f'{pilotpaperpath}/obstable_data.tex', 'w') as fh:
    for key in keys_sorted:
        row = obstable[key]
        fh.write(f"{key:10s} & {row['time']:10.1f} & {row['sessions']:20s} & {row['noiselevel']*1000:10.1f} \\\\\n")
