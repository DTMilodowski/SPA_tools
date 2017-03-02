# This set of functions prepares met data for inclusion into DALEC runs
# Includes ERA interim and TRMM
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

# generate daily met drivers from ERA-Interim and TRMM
# currently do not use shift in time zone, so min, max
# mean and cumulative variables represent GMT "days"
def generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start_date, end_date):
    
   # read in ERA interim.
    ERA_start = '2011-01-01'
    ERA_start_date = np.datetime64(ERA_start,'D')# - np.timedelta64(time_zone_shift,'h')

    sw_eraint, tair_eraint, pptn_eraint, rh_eraint, sp_eraint = wg.readInEraInterim(ERA_file,1) 
    sw_era = sw_eraint[:sw_eraint.size/2].astype('float')  
    pptn_era = pptn_eraint[:sw_eraint.size/2].astype('float')
    airT_era = tair_eraint.astype('float')
    rh_era = rh_eraint.astype('float')
    sp_era = sp_eraint.astype('float')
    vpd_era=(0.6108*np.exp((17.27*airT_era)/(237.3+airT_era)))*(1-rh_era/100)

    n_days_eraint = tair_eraint.size/4
    n_tsteps_ERA = int(n_days_eraint*24/tstep)
    
    mn2t = np.zeros(n_days_eraint)
    mx2t = np.zeros(n_days_eraint)
    ssrd = np.zeros(n_days_eraint)
    vpd = np.zeros(n_days_eraint)
    ERA_dates = np.zeros(n_days_eraint,dtype='datetime64[D]')
    iter_2 = 0
    iter_4 = 0
    for dd in range(0,n_days_eraint):
        # airT and rh both reported quarterly
        mn2t[dd]=np.min(airT_era[dd*iter_4:(dd+1)*iter_4])
        mx2t[dd]=np.max(airT_era[dd*iter_4:(dd+1)*iter_4])
        vpd_era[dd] = np.mean(vpd_era[dd*iter_4:(dd+1)*iter_4])

        # sw rad reported as 12hr totals
        ssrd[dd]=np.sum(sw_era[dd*iter_2:(dd+1)*iter_2])
    
        ERA_dates[dd] = ERA_start+np.timedelta64(1,'D')

    ###############
    # Now load TRMM - TRMM is reported as three hourly totals
    print TRMM_file
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)

    TRMM_dates_only =  TRMM_dates_init.astype('datetime64[D]')
    TRMM_dates = np.unique(TRMM_dates_only)
    N_TRMM = TRMM_dates.size
    TRMM_pptn = np.zeros(N_TRMM)
    for dd in range(0,N_TRMM):
        TRMM_pptn[dd] = np.sum(TRMM_pptn_init[TRMM_dates_init==TRMM_dates[dd]])

    # Now cut time series to period of interest
    start = np.datetime64(start_date,'D')
    end = np.datetime64(end,'D')
    if TRMM_dates[0]>start:
        print "Issue with TRMM - time series starts too late for period of interest"
    if TRMM_dates[-1]<end:
        print "Issue with TRMM - time series finishes too early for period of interest"
    if ERA_dates[0]>start:
        print "Issue with ERA-Interim - time series starts too late for period of interest"
    if ERA_dates[-1]<end:
        print "Issue with ERA-Interim - time series finishes too early for period of interest"
    TRMM_mask = np.all((TRMM_dates>=start,TRMM_dates<=end),axis=0)
    ERA_mask = np.all((ERA_dates>=start,ERA_dates<=end),axis=0)
    
    dates_out = ERA_dates[ERA_mask]
    ssrd_out = ssrd[ERA_mask]
    vpd_out = vpd[ERA_mask]
    mn2t_out = mn2t[ERA_mask]
    mx2t_out = mx2t[ERA_mask]
    pptn_out = TRMM_pptn[TRMM_mask]
    return dates_out,mn2t_out,mx2t_out,vpd_out,ssrd_out,pptn_out
