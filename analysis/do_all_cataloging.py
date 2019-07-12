"""
Estimated execution time: XXX minutes on a 2013 Macbook Pro
"""
import runpy
runpy.run_path('dendro_catalog.py', run_name="__main__") # about 30-40m
runpy.run_path('crossmatch.py', run_name="__main__")
runpy.run_path('fit_gaussians.py', run_name="__main__")
runpy.run_path('concatenate_tables.py', run_name="__main__")
runpy.run_path('crossmatch_plots.py', run_name="__main__")
runpy.run_path('hii_calculations.py', run_name="__main__")
runpy.run_path('sedplots_and_overlays.py', run_name="__main__")
runpy.run_path('overview_plots.py', run_name="__main__")
