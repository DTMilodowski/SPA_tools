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
    print "\t\t remote sensing data"
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

    #--------------------------------------------------------------------------------------------
    print "\t6 - transferring data into output dictionaries"
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

#---------------------------------------------------------------------------------------------
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
def locate_metdata_gaps_using_soil_moisture_time_series(met_data, soil_data, minimum_pptn_rate, STA_LTA_threshold):
    gaps = {}

    met_vars = met_data.keys() 
    N_vars = len(met_vars)
    N_tsteps = met_data[met_vars[0]].size
    template = np.zeros(N_tsteps)

    #---------------------------------------------------------------------------------------------
    # first fdeal with type one gaps
    for vv in range(0,N_vars):
        gaps[met_variables[vv]]=template.copy()
        nodata_mask = np.isnan(met_data[met_vars[vv]])
        gaps[met_variables[vv]][nodata_mask]=1
    
    #---------------------------------------------------------------------------------------------
    # now deal with type two gaps
    # - first of all find out how many sensors are recording. rainfall events must
    #   be deteected across all active sensors to limit false positives.
    N_active_sensors = np.sum((~np.isnan(soil_data['soil_moisture_05cm']),~np.isnan(soil_data['soil_moisture_10cm']),~np.isnan(soil_data['soil_moisture_20cm'])),axis=0)

    # - now use moving averages to define LTA and STA
    STA_half_width = 2. # 2 timesteps (1 hrs) either side of target point - moving average over 2.5 hrs
    LTA_half_width = 48. # 48 timesteps (24hrs) either side of target point - moving average over 48.5 hours
    soil1_STA = met.moving_average_1D(soil_data['soil_moisture_05cm'],STA_half_width)
    soil2_STA = met.moving_average_1D(soil_data['soil_moisture_10cm'],STA_half_width)
    soil3_STA = met.moving_average_1D(soil_data['soil_moisture_20cm'],STA_half_width)

    soil1_LTA = met.moving_average_1D(soil_data['soil_moisture_05cm'],LTA_half_width)
    soil2_LTA = met.moving_average_1D(soil_data['soil_moisture_10cm'],LTA_half_width)
    soil3_LTA = met.moving_average_1D(soil_data['soil_moisture_20cm'],LTA_half_width)

    soil1_filt = soil1_STA/soil1_LTA
    soil2_filt = soil2_STA/soil2_LTA
    soil3_filt = soil3_STA/soil3_LTA

    # find pptn events according to STA/LTA - two options: (i) peak detection; <<(ii) periods above threshold>>.
    rain_event_records = np.sum((soil1_filt>STA_LTA_threshold,soil2_filt>STA_LTA_threshold,soil3_filt>STA_LTA_threshold),axis=0)
    rain_event = rain_event_records == N_active_sensors

    # loop through time series - step through daily - and mark days with "missing rainfall" with number 2
    days = met_data['date'].astype('datetime64[D]')
    days_unique = np.unique(days)
    N_days = days_unique.size
    
    rain_detect_flag = 0
    rain_event_flag = 0
    for dd in range(0,N_days):
        day_index = days==days_unique[dd]
        if np.sum(rain_event[day_index]) > 0:
            rain_event_flag = 1
        if np.max(met_data['pptn'][day_index]) > minimum_pptn_rate:
            rain_detect_flag = 1

        if np.all(rain_detect_flag ==1,rain_event_flag = 0):
            if np.min(N_active_sensors[day_index])>0:
                print "!!! ", days_unique[dd], " scheme not working - soil data fails to pick up precipitation event"
        elif np.all(rain_detect_flag == 0, rain_event_flag == 1):
            gaps['pptn'][day_index]==2    

    rain_detect_flag = 0
    rain_event_flag = 0       
    plot_pptn_detection(met_data['date'], met_data['pptn'], soil_data['soil_moisture_05cm'],soil_data['soil_moisture_10cm'],soil_data['soil_moisture_20cm'], soil1_filt, soil2_filt, soil3_filt, gaps['pptn'])
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

def plot_pptn_detection(dates, pptn, soil1, soil2, soil3, soil1_STA_LTA, soil2_STA_LTA, soil3_STA_LTA, gaps):

    rain_on = np.zeros(len(dates))
    rain_on[rain_sample>0.5]=1*np.ceil(np.max(pptn[np.isfinite(pptn)]))
    nodata=np.zeros(len(dates))
    nodata[np.isnan(pptn)]=1*np.ceil(np.max(rain_sample[np.isfinite(rain_sample)]))
    type1_gap = np.zeros(len(dates))
    type2_gap = np.zeros(len(dates))
    type1_gap[gaps==1] = np.ceil(np.max(pptn[np.isfinite(pptn)]))
    type2_gap[gaps==2] = np.ceil(np.max(pptn[np.isfinite(pptn)]))
    zero = np.zeros(len(dates))

    x_range = np.arange(0,len(dates))
    x_range=x_range/48.

    plt.figure(4, facecolor='White',figsize=[15,8])  
    ax1 = plt.subplot2grid((2,1),(0,0))
    ax1.plot(x_range,pptn,'-',color='#1E7581')
    ax1.fill_between(x_range,zero,type1_gap,color='orange', alpha = 0.2)
    ax1.fill_between(x_range,zero,type2_gap,color='red', alpha = 0.2)
    ax1.fill_between(x_range,zero,pptn,color='#1E7581')
    ax1.set_ylabel('precipitation / mm')
    ax1.set_xlabel('Days since ' + start_date_str)
    ax1.set_ylim(ymax=np.ceil(np.max(pptn[np.isfinite(pptn)])))
    ax1.set_xlim(xmin=0,xmax=np.max(x_range))

    ax2 = plt.subplot2grid((2,1),(1,0),sharex=ax1)
    ax2.plot(x_range,soil1,':',color='#6ACB7B')
    ax2.plot(x_range,soil2,':',color='#24A23A')
    ax2.plot(x_range,soil3,':',color='#006612')
    ax2.plot(x_range,soil1_STA_LTA,'-',color='#6ACB7B',label="5 cm")
    ax2.plot(x_range,soil2_STA_LTA,'-',color='#24A23A',label="10 cm")
    ax2.plot(x_range,soil3_STA_LTA,'-',color='#006612',label="20 cm")
    ax2.fill_between(x_range,zero,type1_gap,color='orange', alpha = 0.2)
    ax2.fill_between(x_range,zero,type2_gap,color='red', alpha = 0.2)
    ax2.set_ylabel('Volumetric soil moisture content')
    ax2.set_xlim(xmin=0,xmax=np.max(x_range))
    ax2.set_ylim(ymin=0,ymax=0.7)
    ax2.legend(loc=2)
    plt.tight_layout(pad=2)
    plt.show()

    return 0
