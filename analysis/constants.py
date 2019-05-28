from astropy import units as u
from radio_beam import Beam

mustang_beam_fwhm = 10*u.arcsec
mustang_central_frequency = 91.5*u.GHz

mgps_beam = Beam(mustang_beam_fwhm)
