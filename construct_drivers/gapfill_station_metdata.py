import numpy as np
from matplotlib import pyplot as plt
import datetime as dt
from matplotlib import rcParams
from scipy import signal

import sys
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/gapfilling/")

import weather_generator_BALI as wg
from statistics_tools import *
import metdata_processing as met
import TRMM_processing as TRMM

# function to read in the met data, soil data and remote sensing data, so that they are all hosted
# in equivalent time series that can be cross-interrogated in preparation for gapfilling.
# The timestep is set in hours and by default is 0.5, which is ideal for running SPA
def load_all_metdata(met_file, soil_file, ERA_file, TRMM_file, start_date, end_date, tstep = 0.5):
    meteorological_data_dict = {}
    soil_data_dict = {} 
    RS_data_dict = {}

    #---------------------------------------------------------------------------------------------
    # load local meteorological data
    met_dates, met_data = met.read_atm_metfile(met_file)
    # load local soil data
    soil_dates, soil_data= met.read_soil_metfile(soil_file)

    earliest_day = np.min(np.min([soil_dates,met_dates]))
    latest_day = np.max(np.max([soil_dates,met_dates]))

    #earliest_day = np.min(met_dates)
    #latest_day = np.max(met_dates)

    dates_series, met_data_host = met.convert_atm_metdata_to_timeseries(met_dates,met_data,earliest_day,latest_day)
    temp, soil_data_host = met.convert_atm_metdata_to_timeseries(soil_dates,soil_data,earliest_day,latest_day)
    
    # shift time series for both so that time is in GMT for comparison with RS products
    time_diff = 8
    dates_series = dates_series - np.timedelta64(time_diff, 'h')

    #---------------------------------------------------------------------------------------------
    # Now read in ERA-Interim data
    # read in ERA interim.
    ERA_start = '01/01/2011 00:00'
    ERA_start_date = np.datetime64(dt.datetime.strptime(ERA_start, '%d/%m/%Y %H:%M'))

    sw_eraint, tair_eraint, pptn_eraint, rh_eraint, sp_eraint = wg.readInEraInterim(ERA_file,1) 
    sw_era = sw_eraint[:sw_eraint.size/2].astype('float')  
    pptn_era = pptn_eraint[:sw_eraint.size/2].astype('float')
    airT_era = tair_eraint.astype('float')
    rh_era = rh_eraint.astype('float')
    sp_era = sp_eraint.astype('float')
    vpd_era=(0.6108*np.exp((17.27*airT_era)/(237.3+airT_era)))*(1-rh_era/100)

    n_days_eraint = tair_eraint.size/4
    n_tsteps_ERA = int(n_days_eraint*24/tstep)
    
    ERA_dates = np.zeros(n_tsteps_ERA).astype('datetime64[m]')
    for i in range(0,n_tsteps_ERA):
        time_increment = int(i*(tstep*60))
        ERA_dates[i]=ERA_start_date+np.timedelta64(time_increment,'m')

    
    #---------------------------------------------------------------------------------------------
    # Now read in TRMM
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)
    TRMM_pptn = TRMM.calc_pptn_specify_tstep(TRMM_pptn_init,tstep)

    n_tsteps_TRMM = TRMM_pptn.size
    TRMM_dates = np.zeross(n_tsteps_TRMM).astype('datetime64[m]')
    for i in range(0,n_tsteps_ERA):
        time_increment = int(i*(tstep*60))
        TRMM_dates[i]=TRMM_dates_init[0]+np.timedelta64(time_increment,'m')

    #---------------------------------------------------------------------------------------------
    # Now use average daily climatology to temporally downscale RS data to specified timestep
    
    # i) Get "average" conditions from met station
    time_in_hours, average_day_pptn, average_day_airT, average_day_RH, average_day_VPD, average_day_PAR, = met.retrieve_monthly_climatology(dates_series,met_data_host)

    # ii) filter these average signals using a moving average
    window_half_width = 1
    boundary_flag = 1 # periodic
    annual_average_day_pptn=met.moving_average_1D(np.mean(average_day_pptn,axis=0),window_half_width,boundary_flag)
    annual_average_day_airT=met.moving_average_1D(np.mean(average_day_airT,axis=0),window_half_width,boundary_flag)
    annual_average_day_RH=met.moving_average_1D(np.mean(average_day_RH,axis=0),window_half_width,boundary_flag)
    annual_average_day_PAR=met.moving_average_1D(np.mean(average_day_PAR,axis=0),window_half_width,boundary_flag)
    annual_average_day_vpd=met.moving_average_1D(np.mean(average_day_VPD,axis=0),window_half_width,boundary_flag)

    # iii) use these averages to model RS metdata at specified timestep
    sw_mod = wg.calc_sw_specify_tstep_with_climatology(sw_era, tstep, annual_average_day_PAR)
    PAR_mod = sw_mod*2.3
    airT_mod = wg.calc_airT_specify_tstep_with_climatology(airT_era, tstep, annual_average_day_airT)
    rh_mod = wg.calc_rh_specify_tstep_with_climatology(rh_era, tstep, annual_average_day_RH)
    vpd_mod =wg.calc_airT_specify_tstep_with_climatology(vpd_era, tstep, annual_average_day_vpd)
    sp_mod = wg.calc_sp_specify_tstep_linear(sp_era, tstep)

    #---------------------------------------------------------------------------------------------
    # Now the data has been prepared, position data on equivalent time series clipped to specified
    # period of interest
    start = np.datetime64(dt.datetime.strptime(start_date, '%d/%m/%Y %H:%M'))
    end = np.datetime64(dt.datetime.strptime(end_date, '%d/%m/%Y %H:%M'))
    output_time_series = np.arange(start,end,np,timedelta64(30,'m'))
    N_tsteps_out = output_time_series.size

    # First deal with the local observations
    airT_station = np.zeros(N_tsteps_out)*np.nan
    pptn_station = np.zeros(N_tsteps_out)*np.nan
    rh_station = np.zeros(N_tsteps_out)*np.nan
    PAR_station = np.zeros(N_tsteps_out)*np.nan
    vpd_station = np.zeros(N_tsteps_out)*np.nan
    BP_station = np.zeros(N_tsteps_out)*np.nan
    swr_station = np.zeros(N_tsteps_out)*np.nan

    soil_moisture_05cm = np.zeros(N_tsteps_out)*np.nan
    soil_moisture_10cm = np.zeros(N_tsteps_out)*np.nan
    soil_moisture_20cm = np.zeros(N_tsteps_out)*np.nan
    
    for tt in range(0,N_tsteps_out):
        if output_time_series[tt] in dates_series:
            index = dates_series == output_time_series[tt]
            airT_station[tt] = met_data_host[index,0]
            pptn_station[tt] = met_data_host[index,1]
            rh_station[tt] = met_data_host[index,5]
            PAR_station[tt] = met_data_host[index,3]
            vpd_station[tt] = met_data_host[index,2]
            BP_station[tt] = met_data_host[index,6]
            swr_station[tt] = met_data_host[index,4]

            soil_moisture_05cm[tt] = soil_data_host[index,6]
            soil_moisture_10cm[tt] = soil_data_host[index,7]
            soil_moisture_20cm[tt] = soil_data_host[index,8]

    # Second deal with RS data
    airT_RS = np.zeros(N_tsteps_out)*np.nan
    pptn_RS = np.zeros(N_tsteps_out)*np.nan
    rh_RS = np.zeros(N_tsteps_out)*np.nan
    PAR_RS = np.zeros(N_tsteps_out)*np.nan
    vpd_RS = np.zeros(N_tsteps_out)*np.nan
    BP_RS = np.zeros(N_tsteps_out)*np.nan
    swr_RS = np.zeros(N_tsteps_out)*np.nan

    for tt in range(0,N_tsteps_out):
        if output_time_series[tt] in ERA_dates:
            index = ERA_dates == output_time_series[tt]
            airT_RS[tt] = airT_mod[index]
            rh_RS[tt] = rh_mod[index]
            PAR_RS[tt] = PAR_mod[index]
            vpd_RS[tt] = vpd_mod[index]
            BP_RS[tt] = sp_mod[index]
            swr_RS[tt] = swr_mod[index]
        if output_time_series[tt] in TRMM_dates:
            index = TRMM_dates == output_time_series[tt]
            pptn_RS[tt] = TRMM_pptn[index]

    # now put all arrays into dictionaries so it is easy to access later
    meteorological_data_dict['date'] = output_time_series.copy()
    meteorological_data_dict['airT'] = airT_station.copy()
    meteorological_data_dict['pptn'] = pptn_station.copy()
    meteorological_data_dict['rh'] = rh_station.copy()
    meteorological_data_dict['PAR'] = PAR_station.copy()
    meteorological_data_dict['vpd'] = vpd_station.copy()
    meteorological_data_dict['BP'] = BP_station.copy()
    meteorological_data_dict['swr'] = swr_station.copy()
    
    soil_meteorological_data_dict['date'] = output_time_series.copy()
    soil_data_dict['soil_moisture_05cm'] = soil_moisture_05cm.copy()    
    soil_data_dict['soil_moisture_10cm'] = soil_moisture_10cm.copy()    
    soil_data_dict['soil_moisture_20cm'] = soil_moisture_20cm.copy()
    
    RS_data_dict['date'] = output_time_series.copy()
    RS_data_dict['airT'] = airT_RS.copy()
    RS_data_dict['pptn'] = pptn_RS.copy()
    RS_data_dict['rh'] = rh_RS.copy()
    RS_data_dict['PAR'] = PAR_RS.copy()
    RS_data_dict['vpd'] = vpd_RS.copy()
    RS_data_dict['BP'] = BP_RS.copy()
    RS_data_dict['swr'] = swr_RS.copy()

    return met_data_dict, soil_data_dict, RS_data_dict


# function to search through the local met and soil data, looking for gaps
#    Two types of gap here:
#    (i)  data gaps where instruments were not recording, due to malfunction, power supply issues etc.
#         These gaps are easily identifiable as having nodata.
#    (ii) periods where the rain guage was recording, but was blocked, either due to leaf fall, which
#         is unsurpisingly common in a tropical rainforest, or poo, prevailingly from birds.  These
#         gaps are more difficult to discern because they are recorded as false zeros.  This next
#         section tries to locate these based on the colocated soil moisture data.
#         Rainfall events are detected based on a STA/LTA threshold.  Gaps are located when no rain is received
#         in the same day that there is a trigger event recorded in the soil moisture time series.
# Type one and type two gaps are returned with a marker specifying which type of gap are present
def locate_metdata_gaps_using_soil_moisture_time_series(met_data, soil_data):
    gaps = {}

    met_variables = met_data.keys() 
    N_vars = len(met_variables)
    N_tsteps = met
    for vv in range(0,N_vars):
        gaps[met_variables[vv]]=0

    
    return gaps

# function to gapfill metdata using remote sensed data
def gapfill_metdata(met_data,RS_data,gaps):
    met_variables = met_data.keys() 
    N_vars = len(met_variables)
    for vv in range(0,N_vars):
        gap_mask = gaps[met_variables[vv]]>0
        met_data[met_variables[vv]][mask]=RS_data[met_variables[vv]][mask]
    return met_data


# function to bias-correct met data.
def bias_correction_monthly(met_data,RS_data):
    bias_corrected_data = {}
    return bias_corrected_data
