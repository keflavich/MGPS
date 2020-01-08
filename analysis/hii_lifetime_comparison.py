from astroquery.vizier import Vizier
import regions
from astropy import wcs
from astropy import coordinates as coord
from astropy.io import fits
from astropy.table import Table
from files import files
from paths import catalog_path

Vizier.find_catalogs('Kalcheva')

Vizier.ROW_LIMIT=1e9
cats = Vizier.get_catalogs('J/A+A/615/A103')
cat = cats[0]

coords = coord.SkyCoord(cat['RAJ2000'], cat['DEJ2000'], frame='fk5')

reglist = regions.io.read_ds9('cutout_regions.reg')

full_table = Table.read(f'{catalog_path}/concatenated_catalog.ipac', format='ascii.ipac')

ttl_uchii = 0
ttl_compact = 0

for regname,fn in files.items():

    reg = [x for x in reglist if x.meta['text'] == regname][0]

    match_cornish = reg.contains(coords, wcs.WCS(fits.getheader(fn)))
    ncornish = match_cornish.sum()

    field_id = full_table['FieldID']
    ncand_field = ((field_id == regname) & (full_table['mgps_and_mm_detected_cm_nondetected'] == 'True')).sum()
    ncompactcand_field = ((field_id == regname) & (full_table['mgps_and_mm_detected_nocm_compact'] == 'True')).sum()
    compact_and_cornish = ((field_id == regname) & (full_table['mgps_and_mm_detected_nocm_compact'] == 'True') & (full_table['CORNISH_EM'] > 0)).sum()

    print(regname, ncornish, ncand_field, ncompactcand_field, compact_and_cornish)

    ttl_uchii += ncornish
    ttl_compact += ncompactcand_field

print(f"Total: CORNISH = {ttl_uchii}   compact = {ttl_compact}   frac = {ttl_compact/ttl_uchii}")
