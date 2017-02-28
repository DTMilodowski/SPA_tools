# This code contains functions to set up a CARDAMOM run for an individual GEM plot, taking in estimates of:
# 1) Biomass - from GEM census
# 2) Root stocks - from GEM plots
# 3) Litter traps - from GEM plots
# 4) LAI estimates - from  MODIS time series
# 5) Soil carbon - HWSD in absence of field data
# 6) Soil respiration? Don't have raw data at present

# These all provide priors into the simulation which in turn help to refine the acceptable parameter space for
# the various parameters in DALEC that will form the starting point for subsequent more detailed SPA modelling.

# CARDAMOM-DALEC is driven with remotely sensed meteorologoical data - from ERA-Interim Reanalysis and TRMM.
# At the plot level, I do not generally include fire and/or deforestation unless directly observed in the
# field campaigns - this is not observed at BALI sites

# The code uses rpy2 (or r2py?) to interface with routines written in R by Luke Smallman to run
# CARDAMOM.  This enables the user to process priors and drivers using python scripts and feed these into
# the preexisting R-routines for running the model.  The advantage of this is that it limits the extent to
# which code is duplicated while retaining flexibility to utilise python for data i/o and processing and display
# of inputs/outputs.

# import some libraries - update as required
import numpy as np
import sys
sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/")
import load_field_data as field

sys.path.append("/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/construct_drivers/")
import construct_met_drivers as construct_met

#---------------------------------------------------------------------------------------------------------------
# First of all, define the relevant input files
# MET DATA
#    -ERA Interim
ERA_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/ERAinterim/BALI_Met/BALI_ERA_Interim_Met.txt'
#    -TRMM
 TRMM_file = '/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/TRMM/g4.areaAvgTimeSeries.TRMM_3B42_007_precipitation.20110101-20160429.117E_4N_117E_4N.csv'

# LAI DATA
#    -MODIS


# FIELD DATA
#    -Plot census data
#    -Fine roots data
#    -Litter data
#    -LiDAR data (for LAI)



#---------------------------------------------------------------------------------------------------------------
# Now get some basic parameters for the run
start_date= '01/01/2011 00:00'
end_date= '01/01/2016 00:00'
timestep = 48*30.





#---------------------------------------------------------------------------------------------------------------
# Process met data
#   -this function returns a dictionary with time series of meteorological variables to be assimilated into
#    CARDAMOM
gapfilled_data = construct_met.EO_metdata(met_station_file,ERA_file,TRMM_file,start_date,end_date)

    start_date_for_cycles = '01/01/2012 00:00'
    end_date_for_cycles = '01/01/2016 00:00'
    gapfilled_data_for_cycles = gapfill_metdata(met_station_file,ERA_file,TRMM_file,start_date_for_cycles,end_date_for_cycles)
    gapfilled_data_ext = extend_metdata(gapfilled_data_for_cycles,2)

    # now concatenate met data
    MetData = {}
    MetData['Time']=np.concatenate((gapfilled_data['Time'],gapfilled_data_ext['Time']))
    MetData['airT']=np.concatenate((gapfilled_data['airT'],gapfilled_data_ext['airT']))
    MetData['pptn']=np.concatenate((gapfilled_data['pptn'],gapfilled_data_ext['pptn']))
    MetData['vpd']= np.concatenate((gapfilled_data['vpd'], gapfilled_data_ext['vpd']))
    MetData['par']= np.concatenate((gapfilled_data['par'], gapfilled_data_ext['par']))
    MetData['swr']= np.concatenate((gapfilled_data['swr'], gapfilled_data_ext['swr']))
    MetData['sp']=  np.concatenate((gapfilled_data['sp'],  gapfilled_data_ext['sp']))

    write_metdata_to_SPA_input(SPAMetDir,SPAMetName,MetData,tstep)
