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
import downscale_rs_data_random_forests as rf
# Met station
met_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_AtmMet_data.csv'
soil_file  = '/home/dmilodow/DataStore_DTM/BALI/SAFE_data/SAFE_FluxTower_SoilMet_data.csv'

# ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'

# TRMM
TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

start_date= '01/01/2011 00:00'
end_date= '01/01/2016 00:00'

met_data_dict, soil_data_dict, RS_data_dict = io.load_all_metdata_for_randomforest(met_file, soil_file, ERA_file, TRMM_file, start_date, end_date)

# remove swr and PAR record from station prior to 22/09/2012 as the sensor was behaving oddly
mask = met_data_dict['date']<np.datetime64('2012-09-22 00:00','m')
met_data_dict['PAR'][mask]=np.nan
met_data_dict['swr'][mask]=np.nan 

# find gaps
minimum_pptn_rate = 0.5
STA_LTA_threshold = 1.1
gaps = gap.locate_metdata_gaps_using_soil_moisture_time_series(met_data_dict, soil_data_dict, minimum_pptn_rate, STA_LTA_threshold)

# random forest regression
rs_variables=np.asarray([RS_data_dict['time'],RS_data_dict['day'],RS_data_dict['airT'],RS_data_dict['airT_n'],RS_data_dict['airT_p'],RS_data_dict['pptn'],RS_data_dict['pptn_d'],RS_data_dict['rh'],RS_data_dict['rh_n'],RS_data_dict['rh_p'],RS_data_dict['PAR'],RS_data_dict['vpd'],RS_data_dict['vpd_n'],RS_data_dict['vpd_p'],RS_data_dict['BP'],RS_data_dict['BP_n'],RS_data_dict['BP_p'],RS_data_dict['swr']]).transpose()

rs_rf={}
n_vars = len(met_data_dict.keys())
for i in range(0,len(met_data_dict.keys())):
    var = met_data_dict.keys()[i]
    if var =='date':
        rs_rf[var] = met_data_dict[var].copy()
    else:
        print var
        target_variable=met_data_dict[var].copy()
        target_variable[gaps[var]!=0]=np.nan
        rs_rf[var] = rf.downsample_with_randomforest(rs_variables,target_variable,n_trees_in_forest = 10000, min_samples_leaf = 50, n_cores = 20)

# gapfill
gapfilled_met_data = gap.gapfill_metdata(met_data_dict,rs_rf,gaps)

SPAMetDir = './'#'/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/SPA_input/met/'
SPAMetName_RS = 'BALI_ERAinterim_TRMM_30mins_v2.csv'
SPAMetName_local = 'BALI_gapfilled_met_station_30mins_v2.csv'

gap.write_metdata_to_SPA_input(SPAMetDir,SPAMetName_local,gapfilled_met_data,gaps)

# create list of zeros for gaps in RS data, since these records are continuous.
gaps_RS = {}
gap_vars = gaps.keys()
N = RS_data_dict['date'].size
for vv in range(0,len(gap_vars)):
    gaps_RS[gap_vars[vv]]=np.zeros(N)

gap.write_metdata_to_SPA_input(SPAMetDir,SPAMetName_RS,RS_data_dict,gaps_RS)

