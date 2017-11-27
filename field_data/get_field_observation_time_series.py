import numpy as np
from matplotlib import pyplot as plt

import load_field_data as field

#==================================
# GET PLOT LEVEL TIME SERIES
#----------------------------------
# Get time series of woody biomass
def get_Cwood_ts(census_file,plot):
    
    census = field.collate_plot_level_census_data(census_file)
    
    plot_biomass = np.sum(census[plot]['C_wood'],axis=0)
    collection_date = np.max(census[plot]['CensusDate'],axis=0)

    return collection_date, plot_biomass

# Get time series of fine root biomass
def get_Croot_NPP_ts(roots_file,plot):
    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_collections = rootNPP[plot]['AccumulationDays'].shape[1]
    collection_dates = []
    previous_collection_dates = []
    NPP = []
    for i in range(0,N_collections):
        # check that there is a collection this year
        if np.max(rootNPP[plot]['CollectionDate'][:,i])>np.datetime64('2000-01-01'):
            collection_dates.append(np.max(rootNPP[plot]['CollectionDate'][:,i]))
            previous_collection_dates.append(np.max(rootNPP[plot]['PreviousCollectionDate'][:,i]))
            # avoid nodata for missing/destroyed traps
            jj = np.isfinite(rootNPP[plot]['FineRootNPP'][:,i])
            NPP.append(np.mean(rootNPP[plot]['FineRootNPP'][jj,i]))
        
    return np.asarray(collection_dates), np.asarray(previous_collection_dates), np.asarray(NPP)
    
# Get time series of litter fall
def get_litterfall_ts(litter_file,plot):
    
    litter = field.read_litterfall_data(litter_file)
    collection_dates = []
    previous_collection_dates = []
    litter_fall = []
    for i in range(0,litter[plot]['mLeaves'].shape[1]):
        collection_dates.append(np.max(litter[plot]['CollectionDate'][:,i]))
        previous_collection_dates.append(np.max(litter[plot]['PreviousCollectionDate'][:,i]))
        
        # avoid nodata for missing/destroyed traps
        jj = np.isfinite(litter[plot]['mLeaves'][:,i])
        litter_fall.append(np.mean(litter[plot]['mLeaves'][jj,i]))
        
    return np.asarray(collection_dates), np.asarray(previous_collection_dates), np.asarray(litter_fall)
