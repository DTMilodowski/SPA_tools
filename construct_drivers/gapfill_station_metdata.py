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
from metdata_processing import *
import TRMM_processing as TRMM


def gapfill_metdata(met_file,ERA_file,TRMM_file,start_date,end_date):

    met_dates, met_data = read_atm_metfile(met_file)
    earliest_day = np.min(met_dates)
    latest_day = np.max(met_dates)
    dates_series, met_data_host = convert_atm_metdata_to_timeseries(met_dates,met_data,earliest_day,latest_day)
    time_diff = 8
    dates_series = dates_series - np.timedelta64(time_diff, 'h')
    
    # Get "average" conditions
    time_in_hours, average_day_pptn, average_day_airT, average_day_RH, average_day_VPD, average_day_PAR, = retrieve_monthly_climatology(dates_series,met_data_host)

    # Relevant parameters for ERA and filter with moving average
    window_half_width = 1
    boundary_flag = 1 # periodic
    annual_average_day_pptn=moving_average_1D(np.mean(average_day_pptn,axis=0),window_half_width,boundary_flag)
    annual_average_day_airT=moving_average_1D(np.mean(average_day_airT,axis=0),window_half_width,boundary_flag)
    annual_average_day_RH=moving_average_1D(np.mean(average_day_RH,axis=0),window_half_width,boundary_flag)
    annual_average_day_PAR=moving_average_1D(np.mean(average_day_PAR,axis=0),window_half_width,boundary_flag)
    annual_average_day_vpd=moving_average_1D(np.mean(average_day_VPD,axis=0),window_half_width,boundary_flag)

    # read in ERA interim.
    ERA_start = '01/01/2011 00:00'
    tstep = 0.5
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

    ###############
    # Now load TRMM
    print TRMM_file
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)
    TRMM_pptn = TRMM.calc_pptn_specify_tstep(TRMM_pptn_init,tstep)

    #########################################################
    # MODEL 1 Linear Interpolation of ERA Interim Reanalysis
    #sw_mod = wg.calc_sw_specify_tstep_linear(sw_era, tstep)
    #pptn_mod = wg.calc_pptn_specify_tstep_linear(pptn_era, tstep)
    #airT_mod = wg.calc_airT_specify_tstep_linear(airT_era, tstep)
    #rh_mod = wg.calc_rh_specify_tstep_linear(rh_era, tstep)
    #PAR_mod = sw_mod1*2.3
    #vpd_mod = wg.calc_airT_specify_tstep_linear(vpd_era, tstep)

    ##################################################################
    # MODEL 3 Climatology-based downsampling of ERA Interim Reanalysis
    sw_mod = wg.calc_sw_specify_tstep_with_climatology(sw_era, tstep, annual_average_day_PAR)
    PAR_mod = sw_mod*2.3
    #pptn_mod = wg.calc_pptn_specify_tstep_with_climatology(pptn_era, tstep, annual_average_day_pptn)
    airT_mod = wg.calc_airT_specify_tstep_with_climatology(airT_era, tstep, annual_average_day_airT)
    rh_mod = wg.calc_rh_specify_tstep_with_climatology(rh_era, tstep, annual_average_day_RH)
    vpd_mod =wg.calc_airT_specify_tstep_with_climatology(vpd_era, tstep, annual_average_day_vpd)
    sp_mod = wg.calc_sp_specify_tstep_linear(sp_era, tstep)
    ####################################################################################################
    # Clip data to time period
    sample_start = np.datetime64(dt.datetime.strptime(start_date, '%d/%m/%Y %H:%M'))
    sample_end = np.datetime64(dt.datetime.strptime(end_date, '%d/%m/%Y %H:%M'))

    # First sample from observations
    temp1 = dates_series>=sample_start
    temp2 = dates_series<sample_end
    indices_obs = temp1*temp2
    dates_sample = dates_series[indices_obs]
    airT_sample = met_data_host[indices_obs,0]
    pptn_sample = met_data_host[indices_obs,1]
    rh_sample = met_data_host[indices_obs,5]
    PAR_sample = met_data_host[indices_obs,3]
    vpd_sample = met_data_host[indices_obs,2]
    BP_sample = met_data_host[indices_obs,6]
    swr_sample = met_data_host[indices_obs,4]


    # Now get dates series for ERA-interim data
    ERA_start_str = '01/01/2011 00:00'
    ERA_end_str = '01/03/2016 00:00'       
    ERA_start = np.datetime64(dt.datetime.strptime(ERA_start_str, '%d/%m/%Y %H:%M'))
    ERA_end = np.datetime64(dt.datetime.strptime(ERA_end_str, '%d/%m/%Y %H:%M'))
    ERA_dates_series = np.arange(ERA_start,ERA_end,np.timedelta64(30,'m'))
    temp1 = ERA_dates_series>=sample_start
    temp2 = ERA_dates_series<sample_end
    indices_ERA = temp1*temp2

    PAR_mod_sample=PAR_mod[indices_ERA]
    swr_mod_sample=PAR_mod_sample/2.3    
    #pptn_mod_sample=pptn_mod[indices_ERA]
    airT_mod_sample=airT_mod[indices_ERA]
    rh_mod_sample=rh_mod[indices_ERA]
    vpd_mod_sample=vpd_mod[indices_ERA]
    sp_mod_sample=sp_mod[indices_ERA]

    BP_sample[BP_sample<=96000]=np.nan
    BP_mean = np.mean(BP_sample[np.isfinite(BP_sample)])

    # Now sample TRMM_pptn_record
    TRMM_start = TRMM_dates_init[0]
    TRMM_end = TRMM_dates_init[-1]      
    TRMM_dates = np.arange(TRMM_start,TRMM_end,np.timedelta64(30,'m'))
    temp1 = TRMM_dates>=sample_start
    temp2 = TRMM_dates<sample_end
    indices_TRMM = temp1*temp2
    TRMM_sample=TRMM_pptn[indices_TRMM]

    #########################
    # Gapfill met station
    airT_gapfilled = airT_sample.copy()
    pptn_gapfilled = pptn_sample.copy()
    vpd_gapfilled = vpd_sample.copy()
    par_gapfilled = PAR_sample.copy()
    swr_gapfilled = swr_sample.copy()
    BP_gapfilled = BP_sample.copy()
    """
    airT_gapfilled[np.isnan(airT_gapfilled)] = airT_mod_sample[np.isnan(airT_gapfilled)]
    pptn_gapfilled[np.isnan(pptn_gapfilled)] = TRMM_sample[np.isnan(pptn_gapfilled)]
    vpd_gapfilled[np.isnan(vpd_gapfilled)] = vpd_mod_sample[np.isnan(vpd_gapfilled)]
    par_gapfilled[np.isnan(par_gapfilled)] = PAR_mod_sample[np.isnan(par_gapfilled)]
    swr_gapfilled[np.isnan(swr_gapfilled)] = swr_mod_sample[np.isnan(swr_gapfilled)]
    BP_gapfilled[np.isnan(BP_gapfilled)] = BP_mean
    """

    airT_gapfilled = airT_mod_sample
    pptn_gapfilled = TRMM_sample
    vpd_gapfilled = vpd_mod_sample
    par_gapfilled = PAR_mod_sample
    swr_gapfilled = swr_mod_sample
    BP_gapfilled = sp_mod_sample


    #############################
    # Create output as dictionary
    MetDataOut={}
    MetDataOut['Time']=dates_sample
    MetDataOut['airT']=airT_gapfilled
    #MetDataOut['pptn']=pptn_gapfilled
    MetDataOut['pptn']=TRMM_sample
    MetDataOut['vpd']=vpd_gapfilled
    MetDataOut['par']=par_gapfilled
    MetDataOut['swr']=swr_gapfilled
    MetDataOut['sp']=BP_gapfilled

    return MetDataOut


def EO_metdata(met_file,ERA_file,TRMM_file,start_date,end_date):

    met_dates, met_data = read_atm_metfile(met_file)
    earliest_day = np.min(met_dates)
    latest_day = np.max(met_dates)
    dates_series, met_data_host = convert_atm_metdata_to_timeseries(met_dates,met_data,earliest_day,latest_day)
    time_diff = 8
    dates_series = dates_series - np.timedelta64(time_diff, 'h')
    
    # Get "average" conditions
    time_in_hours, average_day_pptn, average_day_airT, average_day_RH, average_day_VPD, average_day_PAR, = retrieve_monthly_climatology(dates_series,met_data_host)

    # Relevant parameters for ERA and filter with moving average
    window_half_width = 1
    boundary_flag = 1 # periodic
    annual_average_day_pptn=moving_average_1D(np.mean(average_day_pptn,axis=0),window_half_width,boundary_flag)
    annual_average_day_airT=moving_average_1D(np.mean(average_day_airT,axis=0),window_half_width,boundary_flag)
    annual_average_day_RH=moving_average_1D(np.mean(average_day_RH,axis=0),window_half_width,boundary_flag)
    annual_average_day_PAR=moving_average_1D(np.mean(average_day_PAR,axis=0),window_half_width,boundary_flag)
    annual_average_day_vpd=moving_average_1D(np.mean(average_day_VPD,axis=0),window_half_width,boundary_flag)

    # read in ERA interim.
    ERA_start = '01/01/2011 00:00'
    tstep = 0.5
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

    ###############
    # Now load TRMM
    print TRMM_file
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)
    TRMM_pptn = TRMM.calc_pptn_specify_tstep(TRMM_pptn_init,tstep)

    ##################################################################
    # MODEL Climatology-based downsampling of ERA Interim Reanalysis
    sw_mod = wg.calc_sw_specify_tstep_with_climatology(sw_era, tstep, annual_average_day_PAR)
    PAR_mod = sw_mod*2.3
    #pptn_mod = wg.calc_pptn_specify_tstep_with_climatology(pptn_era, tstep, annual_average_day_pptn)
    airT_mod = wg.calc_airT_specify_tstep_with_climatology(airT_era, tstep, annual_average_day_airT)
    rh_mod = wg.calc_rh_specify_tstep_with_climatology(rh_era, tstep, annual_average_day_RH)
    vpd_mod =wg.calc_airT_specify_tstep_with_climatology(vpd_era, tstep, annual_average_day_vpd)
    sp_mod = wg.calc_sp_specify_tstep_linear(sp_era, tstep)
    ####################################################################################################
    # Clip data to time period
    sample_start = np.datetime64(dt.datetime.strptime(start_date, '%d/%m/%Y %H:%M'))
    sample_end = np.datetime64(dt.datetime.strptime(end_date, '%d/%m/%Y %H:%M'))

    # First sample from observations
    temp1 = dates_series>=sample_start
    temp2 = dates_series<sample_end
    indices_obs = temp1*temp2
    dates_sample = dates_series[indices_obs]

    # Now get dates series for ERA-interim data
    ERA_start_str = '01/01/2011 00:00'
    ERA_end_str = '01/03/2016 00:00'       
    ERA_start = np.datetime64(dt.datetime.strptime(ERA_start_str, '%d/%m/%Y %H:%M'))
    ERA_end = np.datetime64(dt.datetime.strptime(ERA_end_str, '%d/%m/%Y %H:%M'))
    ERA_dates_series = np.arange(ERA_start,ERA_end,np.timedelta64(30,'m'))
    temp1 = ERA_dates_series>=sample_start
    temp2 = ERA_dates_series<sample_end
    indices_ERA = temp1*temp2

    # Now sample TRMM_pptn_record
    TRMM_start = TRMM_dates_init[0]
    TRMM_end = TRMM_dates_init[-1]      
    TRMM_dates = np.arange(TRMM_start,TRMM_end,np.timedelta64(30,'m'))
    temp1 = TRMM_dates>=sample_start
    temp2 = TRMM_dates<sample_end
    indices_TRMM = temp1*temp2

    airT = airT_mod[indices_ERA]
    pptn = TRMM_pptn[indices_TRMM]
    vpd = vpd_mod[indices_ERA]
    par = PAR_mod[indices_ERA]
    swr = PAR_mod[indices_ERA]/2.3 
    BP = sp_mod[indices_ERA]


    #############################
    # Create output as dictionary
    MetDataOut={}
    MetDataOut['Time']=dates_sample
    MetDataOut['airT']=airT
    MetDataOut['pptn']=pptn
    MetDataOut['vpd']=vpd
    MetDataOut['par']=par
    MetDataOut['swr']=swr
    MetDataOut['sp']=BP

    return MetDataOut

