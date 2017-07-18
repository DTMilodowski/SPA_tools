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

import gapfill_station_metdata as gap
import metdata_io as io
import downsample_rs_data_randomforests as rf
# Met station
met_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_AtmMet_data.csv'
soil_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_SoilMet_data.csv'

# ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'

# TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

start_date= '01/01/2011 00:00'
end_date= '01/01/2016 00:00'

met_data_dict, soil_data_dict, RS_data_dict = gap.load_all_metdata(met_file, soil_file, ERA_file, TRMM_file, start_date, end_date)

# compile metdata for ingestion into random forest regression model



rs_rf = rf.downsample_with_randomforest(rs_variables,target_variable)

# remove swr and PAR record from station prior to 22/09/2012 as the sensor was behaving oddly
mask = met_data_dict['date']<np.datetime64('2012-09-22 00:00','m')
met_data_dict['PAR'][mask]=np.nan
met_data_dict['swr'][mask]=np.nan 


minimum_pptn_rate = 0.5
STA_LTA_threshold = 1.1
gaps = gap.locate_metdata_gaps_using_soil_moisture_time_series(met_data_dict, soil_data_dict, minimum_pptn_rate, STA_LTA_threshold)



gapfilled_met_data = gap.gapfill_metdata(met_data_dict,RS_data_dict,gaps)

SPAMetDir = './'#'/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/SPA_input/met/'
SPAMetName_RS = 'BALI_ERAinterim_TRMM_30mins_v1.csv'
SPAMetName_local = 'BALI_gapfilled_met_station_30mins_v1.csv'

gap.write_metdata_to_SPA_input(SPAMetDir,SPAMetName_local,gapfilled_met_data,gaps)

# create list of zeros for gaps in RS data, since these records are continuous.
gaps_RS = {}
gap_vars = gaps.keys()
N = RS_data_dict['date'].size
for vv in range(0,len(gap_vars)):
    gaps_RS[gap_vars[vv]]=np.zeros(N)

gap.write_metdata_to_SPA_input(SPAMetDir,SPAMetName_RS,RS_data_dict,gaps_RS)

