# This code sets up the drivers and priors for a CARDAMOM run based on data from GEM plots.
# 1) Biomass - from GEM census
# 2) Root stocks - from GEM plots
# 3) Litter traps - from GEM plots
# 4) LAI estimates - from  MODIS time series

# Currently litter fluxes included as mean fluxes between collection points.  This will be changed in the
# future to account for cumulative fluxes between dates.

# CARDAMOM-DALEC is driven with remotely sensed meteorologoical data - from ERA-Interim Reanalysis and TRMM.
# At the plot level, I do not generally include fire and/or deforestation unless directly observed in the
# field campaigns - this is not observed at BALI sites

# import some libraries - update as required
import numpy as np
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/CARDAMOM/CARDAMOM_setup_drivers/')
import met_setup as met
import MODIS_setup as MODIS
import field_data_setup as field

# First of all, define the relevant input files
# MET DATA
#    -ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'
#    -TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

# LAI DATA
# assume a constant LAI based on LiDAR data

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
end_date= '01/03/2016'
plot = 'Belian'

# Initiate some arrays to host time series
d,m,y = start_date.split('/')
start = np.datetime64(y+'-'+m+'-'+d,'D')
d,m,y = end_date.split('/')
end = np.datetime64(y+'-'+m+'-'+d,'D')
date = np.arange(start,end, dtype = 'datetime64[D]')

N_t = date.size
mn2t_in = np.zeros(N_t)*-9999.
mx2t_in = np.zeros(N_t)*-9999.
vpd_in = np.zeros(N_t)*-9999.
ssrd_in = np.zeros(N_t)*-9999.
pptn_in = np.zeros(N_t)*-9999.

Cwood_in = np.zeros(N_t)*-9999.
Croot_in = np.zeros(N_t)*-9999.
Litter_in = np.zeros(N_t)*-9999.
Litter_std_in = np.zeros(N_t)*-9999.

LAI_in = np.zeros(N_t)*-9999.
LAI_std_in = np.zeros(N_t)*-9999.

#---------------------------------------------------------------------------------------------------------------
# Process met data
#   -this function returns a dictionary with time series of meteorological variables to be assimilated into
#    CARDAMOM
#   - Variable keys: Time, airT, pptn, vpd, par, swr, sp
met_dates, mn2t, mx2t, vpd, ssrd, TRMM_dates, pptn= met.generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start, end)
N_m = met_dates.size
for dd in range(0,N_m):
    mn2t_in[date == met_dates[dd]] = mn2t[dd]
    mx2t_in[date == met_dates[dd]] = mx2t[dd]
    ssrd_in[date == met_dates[dd]] = ssrd[dd]
    vpd_in[date == met_dates[dd]] = vpd[dd]
N_trmm = TRMM_dates.size
for dd in range(0,N_trmm):
    pptn_in[date == TRMM_dates[dd]] = pptn[dd]
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

# Initially assume average flux rates for litter between collection dates 
N_lit=litter_flux.size

for tt in range(0,N_lit):
    indices = np.all((date>=litter_previous_collection_date[tt], date<litter_collection_date[tt]),axis=0)
    n_days = np.float((litter_collection_date[tt]-litter_previous_collection_date[tt])/ np.timedelta64(1, 'D'))
    Litter_in[indices]= litter_flux[tt]/n_days
    Litter_std_in[indices]= litter_std[tt]/n_days
    
    
# write output to file
