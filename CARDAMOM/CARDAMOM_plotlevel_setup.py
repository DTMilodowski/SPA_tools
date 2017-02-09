# This code contains functions to set up a CARDAMOM run for an individual GEM plot, taking in estimates of:
# 1) Biomass - from GEM census
# 2) Root stocks - from GEM plots
# 3) Litter traps - from GEM plots
# 4) LAI estimates - from  MODIS time series
# 5) Soil carbon - HWSD in absence of field data
# 6) Soil respiration? Don't have raw data at present

# These all provide priors into the simulation which in turn help to refine the acceptable parameter space for
# the various parameters in DALEC that will form the starting point for subsequent more detailed SPA modelling.

# CARDAMOM-DALEC is driven with remotely sensed meteorologoical data - from ERA-Interim Reanalysis and TRMM
