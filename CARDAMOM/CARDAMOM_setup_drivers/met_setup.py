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
def generate_daily_met_drivers_ERAinterim_TRMM(ERA_file, TRMM_file, start_date, end_date, time_zone_shift):
    
   # read in ERA interim.
    ERA_start = '01/01/2011 00:00'
    tstep = 0.5
    ERA_start_date = np.datetime64(dt.datetime.strptime(ERA_start, '%d/%m/%Y %H:%M')) - np.timedelta64(time_zone_shift,'h')

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
    ERA_dates = 
    ###############
    # Now load TRMM
    print TRMM_file
    TRMM_dates_init, TRMM_pptn_init = TRMM.read_TRMM_file(TRMM_file)


