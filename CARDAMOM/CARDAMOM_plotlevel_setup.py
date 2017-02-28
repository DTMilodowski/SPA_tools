# This code contains functions to set up a CARDAMOM run for an individual GEM plot, taking in estimates of:
# 1) Biomass - from GEM census
# 2) Root stocks - from GEM plots
# 3) Litter traps - from GEM plots
# 4) LAI estimates - from  MODIS time series
# 5) Soil carbon - HWSD in absence of field data
# 6) Soil respiration? Don't have raw data at present

# These all provide priors into the simulation which in turn help to refine the acceptable parameter space for
# the various parameters in DALEC that will form the starting point for subsequent more detailed SPA modelling.

# CARDAMOM-DALEC is driven with remotely sensed meteorologoical data - from ERA-Interim Reanalysis and TRMM.
# At the plot level, I do not generally include fire and/or deforestation unless directly observed in the
# field campaigns - this is not observed at BALI sites

# The code uses rpy2 (or r2py?) to interface with routines written in R by Luke Smallman to run
# CARDAMOM.  This enables the user to process priors and drivers using python scripts and feed these into
# the preexisting R-routines for running the model.  The advantage of this is that it limits the extent to
# which code is duplicated while retaining flexibility to utilise python for data i/o and processing and display
# of inputs/outputs.

# import some libraries - update as required
import numpy as np
import sys
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/gapfilling/")

import weather_generator_BALI as wg
from statistics_tools import *
from metdata_processing import *
import TRMM_processing as TRMM

sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/")
import load_field_data as field

sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/")
import construct_met_drivers as construct_met

