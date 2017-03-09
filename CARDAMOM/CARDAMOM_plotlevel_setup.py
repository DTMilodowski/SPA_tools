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
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/")
import load_field_data as field

sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/")
import construct_met_drivers as construct_met

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/CARDAMOM/CARDAMOM_setup_drivers/')
import met_setup as met
import MODIS_setup as MODIS
import field_data_setup as field


#---------------------------------------------------------------------------------------------------------------
# First of all, define the relevant input files
# MET DATA
#    -ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'
#    -TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

# LAI DATA
#    -MODIS
MODIS_file = '/home/dmilodow/DataStore_DTM/BALI/MODIS/BALI-GEM-MODIS-MOD15A2H-006-results.csv'

# FIELD DATA
#    -Plot census data
#    -Fine roots data
#    -Litter data
#    -LiDAR data (for LAI)
census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
roots_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_FineRoots_Stock_NPP_RawData.csv'
litter_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_Litterfall_RawData.csv'

#---------------------------------------------------------------------------------------------------------------
# Now get some basic parameters for the run
start_date= '01/01/2011'
end_date= '01/01/2016'
plot = 'Belian'

# Initiate some arrays to host time series
d,m,y = start_date.split('/')
start = np.datetime64(y+'-'+m+'-'+d,'D')
d,m,y = end_date.split('/')
end = np.datetime64(y+'-'+m+'-'+d,'D')
date = np.arange(start,end, dtype = 'datetime64[D]')

N_t = date.size
mn2t_in = np.zeros(N_t)
mx2t_in = np.zeros(N_t)
vpd_in = np.zeros(N_t)
ssrd_in = np.zeros(N_t)
pptn_in = np.zeros(N_t)

Cwood_in = np.zeros(N_t)
Croot_in = np.zeros(N_t)
Litter_in = np.zeros(N_t)

LAI_in = np.zeros(N_t)
LAI_std_in = np.zeros(N_t)

#---------------------------------------------------------------------------------------------------------------
# Process met data
#   -this function returns a dictionary with time series of meteorological variables to be assimilated into
#    CARDAMOM
#   - Variable keys: Time, airT, pptn, vpd, par, swr, sp
met_dates, mn2t, mx2t, vpd, ssrd, pptn= met.generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start, end)
N_m = met_dates.size
for dd in range(0,N_m):
    mn2t_in[date == met_dates[dd]] = mn2t[dd]
    mx2t_in[date == met_dates[dd]] = mx2t[dd]
    ssrd_in[date == met_dates[dd]] = ssrd[dd]
    vpd_in[date == met_dates[dd]] = vpd[dd]
    pptn_in[date == met_dates[dd]] = pptn[dd]
#---------------------------------------------------------------------------------------------------------------
# Process LAI data
#   -this function returns a dictionary with time series of LAI from MODIS to drive DALEC.
LAI_dict = MODIS.load_point_MODIS_LAI_time_series_from_file(MODIS_file)
LAI = LAI_dict['plot']
N_LAI = LAI['date'].size
for dd in range(0,N_LAI):
    LAI_in[date == LAI['date'][dd]] = LAI['LAI'][dd]
    LAI_std_in[date == LAI['date'][dd]] = LAI['LAI_std'][dd]
#---------------------------------------------------------------------------------------------------------------
# Process field data
census_date, Cwood = field.get_Cwood_ts(census_file,plot)
N_c = Cwood.size
for dd in range(0,N_c):
    Cwood_in[date == census_date[dd]] = Cwood[dd]

root_stock_date, Croot, Croot_std = field.get_Croot_ts(roots_file,plot)
N_r = Croot.size
for dd in range(0,N_r):
    Croot_in[date == root_stock_date[dd]] = Croot[dd]


litter_collection_date, litter_previous_collection_date, litter_flux, litter_std = field.get_litterfall_ts(litter_file,plot)
# work out what to do with litter data
