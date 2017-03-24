# This code sets up the drivers and priors for a CARDAMOM run based on data from GEM plots.
# 1) Biomass - from GEM census
# 2) Root stocks - from GEM plots
# 3) Litter traps - from GEM plots
# 4) LAI estimates - from  MODIS time series
# 5) Soil carbon - from pits in GEM plots.

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
plot = ['Belian','LF','B North']
LAI = [6.69,4.78,3.0]
Csoil = [8295.66.0082857,11275.18,3934.03]

# Initiate some arrays to host time series
d,m,y = start_date.split('/')
start = np.datetime64(y+'-'+m+'-'+d,'D')
d,m,y = end_date.split('/')
end = np.datetime64(y+'-'+m+'-'+d,'D')
date = np.arange(start,end, dtype = 'datetime64[D]')

N_t = date.size
mn2t_in = np.zeros(N_t)-9999.
mx2t_in = np.zeros(N_t)-9999.
vpd_in = np.zeros(N_t)-9999.
ssrd_in = np.zeros(N_t)-9999.
pptn_in = np.zeros(N_t)-9999.

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

for pp in range(0,len(plot)):
    print plot[pp]
    Cwood_in = np.zeros(N_t)-9999.
    Croot_in = np.zeros(N_t)-9999.
    Csoil_in = np.zeros(N_t)-9999.
    Litter_in = np.zeros(N_t)-9999.
    Litter_std_in = np.zeros(N_t)-9999.

    #LAI_in = np.zeros(N_t)-9999.
    #LAI_std_in = np.zeros(N_t)-9999.
    
    # soil carbon reported in g/m2
    Csoil_in[0] = Csoil[pp]

    # Cwood reported in kg (for 1ha plot)
    census_date, Cwood = field.get_Cwood_ts(census_file,plot[pp])
    N_c = Cwood.size
    for dd in range(0,N_c):
        Cwood_in[date == census_date[dd]] = Cwood[dd] *1000./10.**4. # convert kg/ha to g/m^2

    # root stocks reported in Mg/ha
    root_stock_date, Croot, Croot_std = field.get_Croot_ts(roots_file,plot[pp])

    N_r = Croot.size
    for dd in range(0,N_r):
        Croot_in[date == root_stock_date[dd]] = Croot[dd]*10**6/10.**4 # convert Mg/ha to g/m^2
    
        # litter fluxes reported in Mg/ha/yr
    litter_collection_date, litter_previous_collection_date, litter_flux, litter_std = field.get_litterfall_ts(litter_file,plot[pp])

    # Initially assume average flux rates for litter between collection dates 
    N_lit=litter_flux.size
    for tt in range(0,N_lit):
        indices = np.all((date>=litter_previous_collection_date[tt], date<litter_collection_date[tt]),axis=0)
        n_days = np.float((litter_collection_date[tt]-litter_previous_collection_date[tt])/ np.timedelta64(1, 'D'))
        Litter_in[indices]= litter_flux[tt]/n_days * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d
        Litter_std_in[indices]= litter_std[tt]/n_days * (10.**6/10.**4/365.25) # convert Mg/ha/yr to g/m2/d

    # Convert nodata to -9999
    Litter_in[np.isnan(Litter_in)]=-9999.
    Litter_std_in[np.isnan(Litter_std_in)]=-9999.
    Croot_in[np.isnan(Croot_in)]=-9999.
    Cwood_in[np.isnan(Cwood_in)]=-9999.
    Csoil_in[np.isnan(Csoil_in)]=-9999.
    
    mx2t_in[np.isnan(mx2t_in)]=-9999.
    mn2t_in[np.isnan(mn2t_in)]=-9999.
    ssrd_in[np.isnan(ssrd_in)]=-9999.
    vpd_in[np.isnan(vpd_in)]=-9999.
    pptn_in[np.isnan(pptn_in)]=-9999.


    # write output to file
    outfile_drivers = "CARDAMOM_met_drivers_"+plot[pp]+".csv"
    out_drivers = open(outfile_drivers,'w')

    outfile_priors = "CARDAMOM_param_priors_"+plot[pp]+".csv"
    out_priors = open(outfile_priors,'w')

    out_drivers.write('timestep_days, date, mn2t, mx2t, vpd, ssrd, pptn, LAI\n')
    out_priors.write('timestep_days, date, Cwood, Croot, Csoil, LitterFlux, LitterFluxStd\n')
    for tt in range(0,N_t):
        out_drivers.write(str(tt) + ',' + str(date[tt]) + ', ' + str(mn2t_in[tt]) + ',' + str(mx2t_in[tt]) + ',' + str(vpd_in[tt]) + ',' + str(ssrd[tt]) + ',' + str(pptn[tt]) + str(LAI[pp])  + '\n')
        out_priors.write(str(tt) + ',' + str(date[tt]) + ', ' + str(Cwood_in[tt]) + ',' + str(Croot_in[tt]) + ',' + str(Csoil_in[tt]) + ',' + str(Litter_in[tt]) + ',' + str(Litter_std_in[tt]) + '\n')
    out_drivers.close()
    out_priors.close()

