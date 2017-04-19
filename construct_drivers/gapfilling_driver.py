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

minimum_pptn_rate = 0.5
STA_LTA_threshold = 4.
gaps = gap.locate_metdata_gaps_using_soil_moisture_time_series(met_data_dict, soil_data_dict, minimum_pptn_rate, STA_LTA_threshold)


gapfilled_met_data_dict = gapfill_metdata(met_data_dict,RS_data_dict,gaps)


