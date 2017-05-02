#This set of functions prepares field data from GEM plots ready for assimilation into CARDAMOM
import numpy as np
import sys
sys.path.append('../field_data/')

import load_field_data as field

from scipy.interpolate import interp1d

# Get time series of woody biomass
def get_Cwood_ts(census_file,plot):
    census = field.collate_plot_level_census_data(census_file)
    plot_biomass = np.sum(census[plot]['C_wood'],axis=0)
    collection_date = np.max(census[plot]['CensusDate'],axis=0)

    return collection_date, plot_biomass

# Get time series of fine root NPP
def get_root_NPP_ts(roots_file,plot):
    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_collections = rootNPP[plot]['AccumulationDays'].shape[1]
    collection_dates = []
    NPP = []
    NPP_std = []
    for i in range(0,N_collections):
        # check that there is a collection this year
        if np.max(rootNPP[plot]['CollectionDate'][:,i])>np.datetime64('2000-01-01'):
            collection_dates.append(np.max(rootNPP[plot]['CollectionDate'][:,i]))
            previous_collection_dates.append(np.max(rootNPP[plot]['PreviousCollectionDate'][:,i]))
            # avoid nodata for missing/destroyed traps
            jj = np.isfinite(rootNPP[plot]['FineRootNPP'][:,i])
            NPP.append(np.mean(rootNPP[plot]['FineRootNPP'][jj,i]))
            NPP_std.append(np.std(rootNPP[plot]['FineRootNPP'][jj,i]))
        
    return np.asarray(collection_dates), np.asarray(previous_collection_dates), np.asarray(NPP), np.asarray(NPP_std)
    

# Get time series of fine root carbon
"""
def get_Croot_ts(roots_file,plot):
    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_collections = rootStocks[plot]['FineRootStocks'].shape[1]
    collection_dates = []
    previous_collection_dates = []
    Croot = []
    Croot_std = []
    for i in range(0,N_collections):
        # check that there is a collection this year
        if np.max(rootStocks[plot]['CollectionDate'][:,i])>np.datetime64('2000-01-01'):
            collection_dates.append(np.max(rootStocks[plot]['CollectionDate'][:,i]))
            # avoid nodata for missing/destroyed traps
            jj = np.isfinite(rootStocks[plot]['FineRootStocks'][:,i])
            Croot.append(np.mean(rootStocks[plot]['FineRootStocks'][jj,i]))
            Croot_std.append(np.std(rootStocks[plot]['FineRootStocks'][jj,i]))
        
    return np.asarray(collection_dates), np.asarray(Croot), np.asarray(Croot_std)
"""
def get_Croot_ts(roots_file,plot):

    rootStocks,rootNPP = field.read_soil_stocks_and_npp(roots_file)
    N_sp,N_collections = rootStocks[plot]['FineRootStocks'].shape
    collection = []
    for i in range(0,N_collections):
        # check that there is a collection this year
        if np.max(rootStocks[plot]['CollectionDate'][:,i])>np.datetime64('2000-01-01'):
            collection.append(np.max(rootStocks[plot]['CollectionDate'][:,i]))

    collection_dates = np.asarray(collection)
    N_dates = collection_dates.size
    interval = np.zeros(N_dates)
    interval[1:] = interval[1:]-interval[:-1]
    
    Croot_gapfilled = np.zeros((N_sp,N_dates))

    for ss in range(0,N_sp):
        mask = np.isfinite(rootStocks[plot]['FineRootStocks'][ss,:])
        interval_nogaps = np.asarray(interval[mask],dtype='float')
        Croot_nogaps = litter[plot]['LAI'][ss,mask]
        
        f = interp1d(interval_nogaps, Croot_nogaps, kind='cubic')  # without specifying "kind", default is linear

        Croot_gapfilled[ss,:] = f(interval)
    
    Croot_plot_ts = np.mean(Croot_gapfilled,axis=0)
    Croot_std = np.std(Croot_gapfilled,axis=0)
    return collection_dates, Croot_plot_ts, Croot_std



"""
# Get time series of litter fall
def get_litterfall_ts(litter_file,plot):
    
    litter = field.read_litterfall_data(litter_file)
    collection_dates = []
    previous_collection_dates = []
    litter_fall = []
    litter_std = []
    for i in range(0,litter[plot]['rTotal'].shape[1]):
        collection_dates.append(np.max(litter[plot]['CollectionDate'][:,i]))
        previous_collection_dates.append(np.max(litter[plot]['PreviousCollectionDate'][:,i]))
        
        # avoid nodata for missing/destroyed traps
        jj = np.isfinite(litter[plot]['rTotal'][:,i])
        if jj.sum()>0:
            litter_fall.append(np.mean(litter[plot]['rTotal'][jj,i]))
            litter_std.append(np.std(litter[plot]['rTotal'][jj,i]))
        else:
            litter_fall.append(np.nan)
            litter_std.append(np.nan)

    return np.asarray(collection_dates), np.asarray(previous_collection_dates), np.asarray(litter_fall), np.asarray(litter_std)
"""
def get_litterfall_ts(litter_file,plot):
    
    litter = field.read_litterfall_data(litter_file)
    N_sp,N_dates = litter[plot]['rTotal'].shape
    
    collection = []
    previous_collection = []
    for i in range(0,N_dates):
        collection.append(np.max(litter[plot]['CollectionDate'][:,i]))
        previous_collection.append(np.max(litter[plot]['PreviousCollectionDate'][:,i]))

    collection_dates = np.asarray(collection)
    previous_collection_dates = np.asarray(previous_collection)
    interval = np.asarray(collection_dates-previous_collection_dates,dtype='float')
    
    litter_gapfilled = np.zeros((N_sp,N_dates))
    for ss in range(0,N_sp):
        mask = np.isfinite(litter[plot]['LAI'][ss,:])
        interval_nogaps = np.asarray(interval[mask],dtype='float')
        litter_nogaps = litter[plot]['LAI'][ss,mask]
        
        f = interp1d(interval_nogaps, litter_nogaps, kind='cubic')  # without specifying "kind", default is linear

        litter_gapfilled[ss,:] = f(interval)
    
    litter_plot_ts = np.mean(litter_gapfilled,axis=0)
    litter_fall_std = np.std(litter_gapfilled,axis=0)
    return collection_dates, previous_collection_dates, litter_fall_ts, litter_fall_std


# Get time series of LAI using spline interpolation to fill the gaps
def get_LAI_ts(LAI_file,plot):
    LAI = field.load_LAI_time_series(LAI_file)
    N_dates, N_sp = LAI[plot]['LAI'].shape
    interval = np.zeros(N_dates,'timedelta64[D]')
    interval[1:] = LAI[plot]['date'][1:]-LAI[plot]['date'][:-1]

    LAI_gapfilled = np.zeros((N_sp,N_dates))
    for ss in range(0,N_sp):
        mask = np.isfinite(LAI[plot]['LAI'][ss,:])
        interval_nogaps = np.asarray(interval[mask],dtype='float')

        LAI_nogaps = LAI[plot]['LAI'][ss,mask]

        f = interp1d(interval_nogaps, LAI_nogaps, kind='cubic')  # without specifying "kind", default is linear
        LAI_gapfilled[ss,:] = f(interval)

    LAI_plot_ts = np.mean(LAI_gapfilled,axis=0)
    return  LAI[plot]['date'], LAI_plot_ts

