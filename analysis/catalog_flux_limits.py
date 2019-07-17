from astropy import units as u

flux_limits = {20*u.cm: 2*u.mJy, # MAGPIS: Helfand+ 2006
               (20*u.cm).to(u.um): 2*u.mJy, # THOR is 3-7 mJy 7-sigma RMS (Beuther+ 2016) #not used
               6*u.cm: 2.5*u.mJy, # MAGPIS: Giveon+ 2005, CORNISH 2 mjy: Hoare+ 2012
               (6*u.cm).to(u.um): 2.5*u.mJy, # MAGPIS: Giveon+ 2005, CORNISH 2 mjy: Hoare+ 2012
               70*u.um: 20*u.mJy, # from Fig 3 of Molinari 2016, ~20 MJy/sr rms w/6" beam -> 20 mJy/beam
               160*u.um: 26*u.mJy, # from Fig 3 of Molinari 2016, ~10 MJy/sr rms w/10" beam -> 26 mJy/beam
               250*u.um: 85*u.mJy, # 10 MJy/sr rms w/18" beam
               350*u.um: 50*u.mJy, # total guess.  Close, though: 5 MJy/sr w/20" beam gives 53 mJy
               500*u.um: 85*u.mJy, # 2 MJy/sr in 4.24e-8 sr
               870*u.um: 70*u.mJy, # Csengeri+ 2014 fig 1
               1100*u.um: 50*u.mJy, # Ginsburg+ 2013, fig 1
               24*u.um: 200*u.mJy, # total guess, unfortunately: hard to get from Gutermuth+ 2015
              }
