# This function essentially compiles together all of the available LAI data for the BALI
# project for ease of comparison and to facilitate their incorporation into the numerical
# modelling framework.  Sources of LAI data (and other relevant proxies:
# (1)  hemispherical photographs, processed with Hemisfer ->  LAI estimates for each plot
#      have been made at multiple points in time, providing a sparse time series of LAI
# (2)  LiDAR, processed either following Stark et al. (2012) or an alternative approach,
#      (Detto et al., 2015, improved by D T Milodowski et al., in prep).
# (3)  MODIS: 250 m resolution MODIS LAI from the latest MODIS LAI product: MOD15A2H
# (4)  Litter flux estimates
# (5)  Canopy volume estimates.
#----------------------------------------------------------------------------------------
import numpy as np
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/CARDAMOM/CARDAMOM_setup_drivers/')
import met_setup as met
import MODIS_setup as MODIS
import field_data_setup as field

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/')
import inventory_based_LAD_profiles as invent

#-------------------------
# input files
# - field data
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos_TimeSeries.csv'
Litter_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/'
Inventory_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/'

# - lidar canopy profiles
LiDAR_MacHorn_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_LAD_profiles_MacHorn_1m.npz'
LiDAR_rad_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_LAD_profiles_RadiativeTransfer_1m.npz'

# - MODIS data
MODIS_file =  = '/home/dmilodow/DataStore_DTM/BALI/'

# Inventory file, allometry file
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
inventory_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
#-------------------------


# List plots to be assessed
plot = ['Belian','LF','B North','B South', 'E', 'Seraya', 'DC1', 'DC2']
N = len(plot)
N_sp = 25
subplot_area = 20.*20.
# Some dictionaries to host the data
hemiphot = {}
LiDAR_MacHorn = {}
LiDAR_rad = {}
litter = {}
inventory = {}

# Load in the data
# MODIS
MODIS_LAI = MODIS.load_point_MODIS_LAI_time_series_from_file(MODIS_file)

# LiDAR
LiDAR1 = np.load(LiDAR_MacHorn_file)
LiDAR2 = np.load(LiDAR_rad_file)

# Inventory
inventory_data = invent.load_crown_survey_data(inventory_file)
a, b, CF, r_sq, p, H, D = invent.retrieve_crown_allometry(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = invent.calculate_allometric_equations_from_survey(inventory_data)

for pp in range(0,N):
    # hemiphotos -> upload subplot level data
    hemiphot[plot[pp]]={}
    hemiphot[plot[pp]]['date'],hemiphot[plot[pp]]['LAI'] = field.get_subplot_LAI_ts(LAI_file,plot[pp])

    # LiDAR data -> retrieve subplot level data
    LiDAR_MacHorn[plot[pp]]['LAI'] = np.sum(LiDAR1[plot[pp]],axis=1)
    LiDAR_rad[plot[pp]]['LAI'] = np.sum(LiDAR2[plot[pp]][:,:,-1],axis=1)

    # Inventory file -> retrieve subplot level data
    # set up array to host inventory profiles
    inventory_LAI = np.zeros(N_sp)
    for ss in range(0,N_sp):
        subplot = ss+1
        mask = np.all((inventory_data['plot']==plot[pp],inventory_data['subplot']==subplot),axis=0)
                      
        Ht,Area,Depth = invent.calculate_crown_dimensions(inventory_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        inventory_LAD, CanopyV = invent.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)
        inventory_LAI[ss] = np.sum(inventory_LAD)
    inventory['volume'] = inventory_LAI.copy()
    

    # Litter file -> upload subplot level data
    litter[plot[pp]]={}
    litter[plot[pp]]['date'],litter[plot[pp]]['flux'] = field.get_subplot_litterfall_ts(Litter_file,plot[pp])
