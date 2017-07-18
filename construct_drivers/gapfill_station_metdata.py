import numpy as np
from matplotlib import pyplot as plt
import datetime as dt
from matplotlib import rcParams
from scipy import signal
from sklearn.ensemble import RandomForestRegressor

import sys
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/")
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/gapfilling/")

import weather_generator_BALI as wg
from statistics_tools import *
import metdata_processing as met
import TRMM_processing as TRMM


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

    # set up blank gap arrays
    for vv in range(0,N_vars):
        if met_vars[vv]!='date':
            gaps[met_vars[vv]]=template.copy()
    #---------------------------------------------------------------------------------------------
    # deal with type two gaps
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

    peaks1 = find_maxima(soil1_filt, STA_LTA_threshold)
    peaks2 = find_maxima(soil2_filt, STA_LTA_threshold)
    peaks3 = find_maxima(soil3_filt, STA_LTA_threshold)

    # find pptn events according to STA/LTA - two options: (i) peak detection; <<(ii) periods above threshold>>.
    #rain_event_records = np.sum((soil1_filt>STA_LTA_threshold,soil2_filt>STA_LTA_threshold,soil3_filt>STA_LTA_threshold),axis=0)
    rain_event_records = np.sum((peaks1,peaks2,peaks3),axis=0)
    rain_event = np.all((rain_event_records>0,rain_event_records >= N_active_sensors-1,rain_event_records > 0),axis=0)
    #rain_event = rain_event_records>0
    
    # loop through time series - step through daily - and mark days with "missing rainfall" with number 2
    days = met_data['date'].astype('datetime64[D]')
    days_unique = np.unique(days)
    N_days = days_unique.size
    
    rain_detect_flag = False
    rain_event_flag = False
    type_2_gap_continue = False

    for dd in range(0,N_days):
        day_index = days==days_unique[dd]
        if np.sum(rain_event[day_index]) > 0:
            rain_event_flag = True # rain indicated by soil moisture records

        if np.max(met_data['pptn'][day_index]) > minimum_pptn_rate:
            print dd, " halleluia it rained"
            rain_detect_flag = True # rain detected in gauge
            type_2_gap_continue = False # as rain detected, therefore gauge is working again!
            gaps['pptn'][day_index]=0
        elif type_2_gap_continue:
            print dd, " continue type 2 gap"
            gaps['pptn'][day_index]=2    
        elif np.all((~rain_detect_flag, rain_event_flag)):
            print dd, " found type 2 gap"
            gaps['pptn'][day_index]=2    
            type_2_gap_continue = True        
        else:
            print dd, " no rain, but working ok"
            gaps['pptn'][day_index]=0

        """
        if np.all((rain_detect_flag,~rain_event_flag)):
            if np.min(N_active_sensors[day_index])>0:
                print dd,  "!!! ", days_unique[dd], " scheme not working - soil data fails to pick up precipitation event"
        """

        rain_detect_flag = False
        rain_event_flag = False  


    #---------------------------------------------------------------------------------------------
    # now quickly deal with type one gaps - where data was not recorded at all
    for vv in range(0,N_vars):
        if met_vars[vv]!='date':
            nodata_mask = np.isnan(met_data[met_vars[vv]])
            gaps[met_vars[vv]][nodata_mask]=1

    #plot_pptn_detection(met_data['date'], met_data['pptn'], soil_data['soil_moisture_05cm'],soil_data['soil_moisture_10cm'],soil_data['soil_moisture_20cm'], soil1_filt, soil2_filt, soil3_filt, gaps['pptn'])
    
    return gaps

# Simple function to find local maxima based on immediate neighbourhood
def find_maxima(signal_y, threshold=0):
    N=signal_y.size
    peak_present = np.zeros(N)
    for i in range(1,N-1):
        if np.all((signal_y[i]>=threshold,signal_y[i]>signal_y[i-1],signal_y[i]>signal_y[i+1])):
            peak_present[i] = 1
    return peak_present

# function to gapfill metdata using remote sensed data
def gapfill_metdata(met_data,RS_data,gaps):
    met_variables = met_data.keys() 
    N_vars = len(met_variables)
    for vv in range(0,N_vars):
        if met_variables[vv]!='date':
            gap_mask = gaps[met_variables[vv]]>0
            met_data[met_variables[vv]][gap_mask]=RS_data[met_variables[vv]][gap_mask]
    return met_data


# function to bias-correct met data.
def bias_correction_monthly(met_data,RS_data):
    bias_corrected_data = {}
    return bias_corrected_data

def plot_pptn_detection(dates, pptn, soil1, soil2, soil3, soil1_STA_LTA, soil2_STA_LTA, soil3_STA_LTA, gaps):

    import matplotlib.transforms as mtransforms


    rain_on = np.zeros(len(dates))
    rain_on[pptn>0.5]=1*np.ceil(np.max(pptn[np.isfinite(pptn)]))
    type1_gap = np.zeros(len(dates))
    type2_gap = np.zeros(len(dates))
    type1_gap[gaps==1] = np.ceil(np.max(pptn[np.isfinite(pptn)]))
    type2_gap[gaps==2] = np.ceil(np.max(pptn[np.isfinite(pptn)]))
    zero = np.zeros(len(dates))

    x_range = np.arange(0,len(dates))
    x_range=x_range/48.

    plt.figure(4, facecolor='White',figsize=[15,10])  
    ax1 = plt.subplot2grid((3,1),(0,0))
    trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax1.plot(x_range,pptn,'-',color='#1E7581')
    ax1.fill_between(x_range, 0, 1, where=gaps == 2, facecolor='green', alpha=0.5,edgecolor='None', transform=trans)
    ax1.fill_between(x_range, 0, 1, where=gaps == 1, facecolor='red', alpha=0.5, edgecolor='None', transform=trans)
    ax1.fill_between(x_range,zero,pptn,color='#1E7581')
    ax1.set_ylabel('precipitation / mm')
    #ax1.set_xlabel('Days since ' + start_date_str)
    ax1.set_ylim(ymin=0,ymax=np.ceil(np.max(pptn[np.isfinite(pptn)])))
    ax1.set_xlim(xmin=0,xmax=np.max(x_range))

    ax2 = plt.subplot2grid((3,1),(1,0),sharex=ax1)
    trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)
    ax2.plot(x_range,soil1,'-',color='#6ACB7B')
    ax2.plot(x_range,soil2,'-',color='#24A23A')
    ax2.plot(x_range,soil3,'-',color='#006612')
    ax2.fill_between(x_range, 0, 1, where=gaps == 2, facecolor='green', alpha=0.5,edgecolor='None', transform=trans)
    ax2.fill_between(x_range, 0, 1, where=gaps == 1, facecolor='red', alpha=0.5, edgecolor='None', transform=trans)
    ax2.set_ylabel('Volumetric soil moisture content')
    ax2.set_xlim(xmin=0,xmax=np.max(x_range))
    ax2.set_ylim(ymin=0,ymax=0.7)
    ax2.legend(loc=2)

    ax3 = plt.subplot2grid((3,1),(2,0),sharex=ax1)
    trans = mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
    ax3.plot(x_range,soil1_STA_LTA,'-',color='#6ACB7B',label="5 cm")
    ax3.plot(x_range,soil2_STA_LTA,'-',color='#24A23A',label="10 cm")
    ax3.plot(x_range,soil3_STA_LTA,'-',color='#006612',label="20 cm")
    ax3.fill_between(x_range, 0, 1, where=gaps == 2, facecolor='green', alpha=0.5, edgecolor='None',transform=trans)
    ax3.fill_between(x_range, 0, 1, where=gaps == 1, facecolor='red', alpha=0.5,edgecolor='None', transform=trans)
    ax3.set_ylabel('STA/LTA')
    ax3.set_xlim(xmin=0,xmax=np.max(x_range))
    ax3.set_ylim(ymin=np.nanmin((soil1_STA_LTA, soil2_STA_LTA, soil3_STA_LTA)),ymax=np.nanmax((soil1_STA_LTA, soil2_STA_LTA, soil3_STA_LTA)))
    plt.tight_layout(pad=2)
    plt.show()

    return 0
