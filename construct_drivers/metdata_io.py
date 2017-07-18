import numpy as np
import datetime as dt
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
    met_data_dict = {}
    soil_data_dict = {} 
    RS_data_dict = {}
    print "Loading met data"
    #---------------------------------------------------------------------------------------------
    # load local meteorological data
    print "\t1 - loading station data" 
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
    print "\t2 - loading ERA-Interim Reanalyses"
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
    print "\t3 - loading TRMM precipitation data"
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)
    TRMM_pptn = TRMM.calc_pptn_specify_tstep(TRMM_pptn_init,tstep)

    n_tsteps_TRMM = TRMM_pptn.size
    TRMM_dates = np.zeros(n_tsteps_TRMM).astype('datetime64[m]')
    for i in range(0,n_tsteps_TRMM):
        time_increment = int(i*(tstep*60))
        TRMM_dates[i]=TRMM_dates_init[0]+np.timedelta64(time_increment,'m')

    #---------------------------------------------------------------------------------------------
    # Now use average daily climatology to temporally downscale RS data to specified timestep
    print "\t4 - use average day climatology to downsample met data to specified timestep"
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
    swr_mod = wg.calc_sw_specify_tstep_with_climatology(sw_era, tstep, annual_average_day_PAR)
    PAR_mod = swr_mod*2.3
    airT_mod = wg.calc_airT_specify_tstep_with_climatology(airT_era, tstep, annual_average_day_airT)
    rh_mod = wg.calc_rh_specify_tstep_with_climatology(rh_era, tstep, annual_average_day_RH)
    vpd_mod =wg.calc_airT_specify_tstep_with_climatology(vpd_era, tstep, annual_average_day_vpd)
    sp_mod = wg.calc_sp_specify_tstep_linear(sp_era, tstep)

    #---------------------------------------------------------------------------------------------
    # Now the data has been prepared, position data on equivalent time series clipped to specified
    # period of interest
    print "\t5 - clipping time series"
    start = np.datetime64(dt.datetime.strptime(start_date, '%d/%m/%Y %H:%M'))
    end = np.datetime64(dt.datetime.strptime(end_date, '%d/%m/%Y %H:%M'))
    output_time_series = np.arange(start,end,np.timedelta64(30,'m'))
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

    print "\t\t local_observations"
    out_indices = np.all((output_time_series >= dates_series[0], output_time_series <= dates_series[-1]),axis=0)
    data_i1 = 0
    data_i2 = dates_series.size
    data_indices =  np.arange(0,dates_series.size)
    if dates_series[0]<output_time_series[0]:
        data_i1 = data_indices[dates_series==output_time_series[0]][0]
    if dates_series[-1]>output_time_series[-1]:
        data_i2 = data_indices[dates_series==output_time_series[-1]][0]+1

    airT_station[out_indices] = met_data_host[data_i1:data_i2,0]
    pptn_station[out_indices] = met_data_host[data_i1:data_i2,1]
    rh_station[out_indices] = met_data_host[data_i1:data_i2,5]
    PAR_station[out_indices] = met_data_host[data_i1:data_i2,3]
    vpd_station[out_indices] = met_data_host[data_i1:data_i2,2]
    BP_station[out_indices] = met_data_host[data_i1:data_i2,6]
    swr_station[out_indices] = met_data_host[data_i1:data_i2,4]

    soil_moisture_05cm[out_indices] = soil_data_host[data_i1:data_i2,6]
    soil_moisture_10cm[out_indices] = soil_data_host[data_i1:data_i2,7]
    soil_moisture_20cm[out_indices] = soil_data_host[data_i1:data_i2,8]

    # Second deal with RS data
    airT_RS = np.zeros(N_tsteps_out)*np.nan
    pptn_RS = np.zeros(N_tsteps_out)*np.nan
    rh_RS = np.zeros(N_tsteps_out)*np.nan
    PAR_RS = np.zeros(N_tsteps_out)*np.nan
    vpd_RS = np.zeros(N_tsteps_out)*np.nan
    BP_RS = np.zeros(N_tsteps_out)*np.nan
    swr_RS = np.zeros(N_tsteps_out)*np.nan

    print "\t\t remote sensing data"
    out_indices = np.all((output_time_series >= ERA_dates[0], output_time_series <= ERA_dates[-1]),axis=0)
    data_i1 = 0
    data_i2 = ERA_dates.size
    data_indices =  np.arange(0,ERA_dates.size)
    if ERA_dates[0]<output_time_series[0]:
        data_i1 = data_indices[ERA_dates==output_time_series[0]][0]
    if ERA_dates[-1]>output_time_series[-1]:
        data_i2 = data_indices[ERA_dates==output_time_series[-1]][0]+1

    airT_RS[out_indices] = airT_mod[data_i1:data_i2]
    rh_RS[out_indices] = rh_mod[data_i1:data_i2]
    PAR_RS[out_indices] = PAR_mod[data_i1:data_i2]
    vpd_RS[out_indices] = vpd_mod[data_i1:data_i2]
    BP_RS[out_indices] = sp_mod[data_i1:data_i2]
    swr_RS[out_indices] = swr_mod[data_i1:data_i2]

    out_indices = np.all((output_time_series >= TRMM_dates[0], output_time_series <= TRMM_dates[-1]),axis=0)
    data_i1 = 0
    data_i2 = TRMM_dates.size
    data_indices =  np.arange(0,TRMM_dates.size)
    if TRMM_dates[0]<output_time_series[0]:
        data_i1 = data_indices[TRMM_dates==output_time_series[0]][0]
    if TRMM_dates[-1]>output_time_series[-1]:
        data_i2 = data_indices[TRMM_dates==output_time_series[-1]][0]+1
    pptn_RS[out_indices] = TRMM_pptn[data_i1:data_i2]

    #--------------------------------------------------------------------------------------------
    print "\t6 - transferring data into output dictionaries"
    # now put all arrays into dictionaries so it is easy to access later
    met_data_dict['date'] = output_time_series.copy()
    met_data_dict['airT'] = airT_station.copy()
    met_data_dict['pptn'] = pptn_station.copy()
    met_data_dict['rh'] = rh_station.copy()
    met_data_dict['PAR'] = PAR_station.copy()
    met_data_dict['vpd'] = vpd_station.copy()
    met_data_dict['BP'] = BP_station.copy()
    met_data_dict['swr'] = swr_station.copy()
    
    soil_data_dict['date'] = output_time_series.copy()
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


#==========================================================================================
# stuff to prep data for ingestion into random forest algorithm - i.e. daily uniform
# distributions of EO parameters
#------------------------------------
# ERA-Interim datasets
# Temperature and relative humidity and surface pressure
# These are quarterly instantaneous measurements.  For each timestep, useful to know current
# temperature and temperature either side.
# 
def calc_airT_specify_tstep_uniform(airt, tstep):
      n_days_eraint = airt.size/4

      # timesteps at 00:00, 06:00, 12:00, 18:00
      # Therefore, for final timestep, assume that midnight temperature matches previous day
      # and for first timetep, previous temperature matches 18:00 temperature from same day
      airt_prev = np.roll(airt,1)
      airt_prev[0] = airt[3]
      airt_next = np.roll(airt,-1)
      airt_next[-1] = airt[-4]

      n_tsteps = int(n_days_eraint*24/tstep)
      airt_out = np.zeros(n_tsteps)
      airt_p_out = np.zeros(n_tsteps)
      airt_n_out = np.zeros(n_tsteps)

      counter = 0
      index = 0
      tsteps_in_quarter = int(6/tstep)
      
      for tt in range(0,n_tsteps):
            airt_out[tt]=airt[index]
            airt_p_out[tt]=airt_prev[index]
            airt_n_out[tt]=airt_next[index]
            counter+=1
            if counter >= tsteps_in_quarter:
                  # reached end of time segment, so restart counter and increment index
                  index+=1
                  counter=0

      return airt_out, airt_p_out, airt_n_out


def calc_rh_specify_tstep_uniform(rh, tstep):
      n_days_eraint = rh.size/4
      # timesteps at 00:00, 06:00, 12:00, 18:00
      # Therefore, for final timestep, assume that midnight Rh matches previous day
      # and for first timetep, previous Rh matches 18:00 Rh from same day
      rh_prev = np.roll(rh,1)
      rh_prev[0] = rh[3]
      rh_next = np.roll(rh,-1)
      rh_next[-1] = rh[-4]

      n_tsteps = int(n_days_eraint*24/tstep)
      rh_out = np.zeros(n_tsteps)
      rh_p_out = np.zeros(n_tsteps)
      rh_n_out = np.zeros(n_tsteps)

      counter = 0
      index = 0
      tsteps_in_quarter = int(6/tstep)
      
      for tt in range(0,n_tsteps):
            rh_out[tt]=rh[index]
            rh_p_out[tt]=rh_prev[index]
            rh_n_out[tt]=rh_next[index]
            counter+=1
            if counter >= tsteps_in_quarter:
                  # reached end of time segment, so restart counter and increment index
                  index+=1
                  counter=0

      return rh_out,rh_p_out,rh_n_out

def calc_sp_specify_tstep_uniform(sp, tstep):
      n_days_eraint = sp.size/4
      # timesteps at 00:00, 06:00, 12:00, 18:00
      sp_prev = np.roll(sp,1)
      sp_prev[0] = sp[3]
      sp_next = np.roll(sp,-1)
      sp_next[-1] = sp[-4]

      n_tsteps = int(n_days_eraint*24/tstep)
      sp_out = np.zeros(n_tsteps)
      sp_p_out = np.zeros(n_tsteps)
      sp_n_out = np.zeros(n_tsteps)

      counter = 0
      index = 0
      tsteps_in_quarter = int(6/tstep)
      
      for tt in range(0,n_tsteps):
            sp_out[tt]=sp_temp[index]
            sp_p_out[tt]=sp_prev[index]
            sp_n_out[tt]=sp_next[index]
            counter+=1
            if counter >= tsteps_in_quarter:
                  # reached end of time segment, so restart counter and increment index
                  index+=1
                  counter=0

      return sp_out,sp_p_out,sp_n_out

# shortwave radiation
# These are 12 hour cumulative totals, distributed uniformly for ingestion into random forest.
def calc_sw_specify_tstep_uniform(sw, tstep):
      n_days_eraint = sw.size/2

      n_tsteps = int(n_days_eraint*24/tstep)
      sw_out = np.zeros(n_tsteps)

      counter = 0
      index = 0
      tsteps_in_half = int(12/tstep)
      seconds_in_tstep = tstep*60.*60.
      
      for tt in range(0,n_tsteps):
            if counter < tsteps_in_half:
                  sw_out[tt]=sw[index]/(float(tsteps_in_half)*seconds_in_tstep)
                  counter+=1
            else:
                  # reached end of time segment, so restart counter and increment index:
                  counter=0
                  index+=1
                  # divide by number of timesteps in half multiplied by number of seconds in a timestep to give average flux in time segment
                  sw_out[tt]=sw[index]/(float(tsteps_in_half)*seconds_in_tstep)
                  counter+=1
      
      return sw_out

#------------------------------------
# TRMM
# Precipitation likely to be weakly correlated at sub-daily timestep, aggregate
# over 24 hour timestep - precipitation will then be allocated based on other
# data
# input: 3 hourly precipitation (in mm/hr)
# output:3 hourly precipitation  at tstep required, and equivalent but for
# precipitation distributed evenly through each day - I grab both as if one is
# not important, it won't be a strong component of the random forest regressor
def calc_TRMM_specify_tstep_uniform(sw, tstep):
   # convert 3 hourly precipitation (in mm/hr) to specified timestep, output units are mm
    n_tsteps_in = pptn_in.size
    tstep_in = 3.
    # input data at three hours temporal resolution
    n_tsteps_out = int(n_tsteps_in*tstep_in/tstep_out)
    pptn_out = np.zeros(n_tsteps_out)
    pptn_day_out = np.zeros(n_tsteps_out)
    tstep_d = 24/tstep_out
    
    counter = 0.
    index = 0
      
    for tt in range(0,n_tsteps_out):
        if counter < tstep_in:
            pptn_out[tt]=pptn_in[index]*tstep_out
            counter+=tstep_out
        else:
            # reached end of time segment, so restart counter and increment index:
            counter=0.
            index+=1
            pptn_out[tt]=pptn_in[index]*tstep_out
            counter+=tstep_out

    for tt in range(0,n_tsteps_out/tstep_d):
        pptn_day_out[tt*tstep_d:(tt+1)*tstep_d]=np.mean(pptn[tt*tstep_d:(tt+1)*tstep_d])

    return pptn_out, pptn_day_out


#======================================================================================
## DATA OUTPUT
# Write some output files
def write_metdata_to_SPA_input(MetDir,MetName,MetData,gaps,tstep_mins=30.):

    out = open(MetDir+MetName,'w')
    out.write('time_days, date, airt,wind_spd,sw_rad,vpd,ppfd,precip,sfc_pressure, gap_airt, gap_sw, gap_vpd, gap_ppfd, gap_pptn, gap_sfc_pressure\n')

    N_timesteps = MetData['date'].size
    tstep_days = tstep_mins/(60.*24.)
    t=0
    for i in range(0,N_timesteps):
        t += tstep_days
        out.write(str(t) + ',' + str(MetData['date'][i]) + ',' + str(MetData['airT'][i]) + ', 0.2,' + str(MetData['swr'][i]) + ',' + str(MetData['vpd'][i]) + ',' + str(MetData['PAR'][i]) + ',' + str(MetData['pptn'][i]) + ',' + str(MetData['BP'][i]) + ',' + str(gaps['airT'][i]) + ',' + str(gaps['swr'][i]) +  ',' + str(gaps['vpd'][i]) + ',' + str(gaps['PAR'][i]) + ',' + str(gaps['pptn'][i]) + ',' + str(gaps['BP'][i]) +'\n')

    out.close()
    return 0


# Write some output files
def write_metdata_to_CARDAMOM_input(MetDir,MetName,MetData):

    out = open(MetDir+MetName,'w')
    out.write('time_days, date, airt,wind_spd,sw_rad,vpd,ppfd,precip,sfc_pressure, gap_airt, gap_sw, gap_vpd, gap_ppfd, gap_pptn, gap_sfc_pressure\n')

    N_timesteps = MetData['date'].size
    tstep_days = tstep_mins/(60.*24.)
    t=0
    for i in range(0,N_timesteps):
        t += tstep_days
        out.write(str(t) + ',' + str(MetData['date'][i]) + ',' + str(MetData['airT'][i]) + ', 0.2,' + str(MetData['swr'][i]) + ',' + str(MetData['vpd'][i]) + ',' + str(MetData['PAR'][i]) + ',' + str(MetData['pptn'][i]) + ',' + str(MetData['BP'][i]) + ',' + str(gaps['airT'][i]) + ',' + str(gaps['swr'][i]) +  ',' + str(gaps['vpd'][i]) + ',' + str(gaps['PAR'][i]) + ',' + str(gaps['pptn'][i]) + ',' + str(gaps['BP'][i]) +'\n')

    out.close()
    return 0
