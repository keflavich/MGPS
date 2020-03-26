from astropy import units as u
from astropy.io import fits
import files

# these are straight averages of the mean values in Table 3
hdr_data = {'G31': {'BMAJ': 10.0*u.arcsec, 'BMIN':9.1*u.arcsec, 'BEAMAREA': 126*u.arcsec**2,},
            'G12': {'BMAJ': 10.0*u.arcsec, 'BMIN':9.0*u.arcsec, 'BEAMAREA': 126*u.arcsec**2,},
            'G43': {'BMAJ': 9.85*u.arcsec, 'BMIN':9.15*u.arcsec, 'BEAMAREA': 121.5*u.arcsec**2,},
            'G49': {'BMAJ': 9.7*u.arcsec, 'BMIN':9.1*u.arcsec, 'BEAMAREA': 117*u.arcsec**2,},
            'G01': {'BMAJ': 10.3*u.arcsec, 'BMIN':9.1*u.arcsec, 'BEAMAREA': 128*u.arcsec**2,},
            'G29': {'BMAJ': 10.6*u.arcsec, 'BMIN':9.1*u.arcsec, 'BEAMAREA': 130*u.arcsec**2,},
            'G34': {'BMAJ': 10.0*u.arcsec, 'BMIN':9.3*u.arcsec, 'BEAMAREA': 133*u.arcsec**2,},
           }



for key in hdr_data:
    hdr_data[key]['REFFREQ'] = (90*u.GHz).to(u.Hz).value
    hdr_data[key]['JYTOK'] = (1*u.Jy).to(u.K,
                                         u.brightness_temperature(90*u.GHz,
                                                                  hdr_data[key]['BEAMAREA'])).value
    for kk in ['BMAJ','BMIN']:
        hdr_data[key][kk] = hdr_data[key][kk].to(u.deg).value
    hdr_data[key]['BEAMAREA'] = (hdr_data[key]['BEAMAREA'].value, "arcsec^2")

for fkey in files.files:
    fn = files.files[fkey]
    print(fkey,fn)
    fh = fits.open(fn, mode='update')
    for key, value in hdr_data[fkey].items():
        fh[0].header[key] = value
    del fh[0].header['BPA']
    fh.flush()
