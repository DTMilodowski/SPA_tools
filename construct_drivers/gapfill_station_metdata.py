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

# function to read in the met data, soil data and remote sensing data, so that they are all hosted in equivalent time series that can be cross-interrogated in preparation for gapfilling.
def load_all_metdata(met_file, soil_file, ERA_file, TRMM_file):
    met_data = {} 
    soil_data = {} 
    RS_data = {}


    return met_data, soil_data, RS_data


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

    return gaps


# function to bias-correct met data.
def bias_correction_monthly(met_data,RS_data):
    bias_corrected_data = {}
    return bias_corrected_data



# function to gapfill metdata using remote sensed data
def gapfill_metdata(met_data,RS_data,gaps):
    gapfilled_data = {}
    return gapfilled_data
