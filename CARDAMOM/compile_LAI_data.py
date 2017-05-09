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

from matplotlib import pyplot as plt
from matplotlib import rcParams
from datetime import datetime
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

sys.path.append('/home/dmilodow/DataStore_DTM/USEFUL_PYTHON_STUFF/generic_plotting/src/')
import violin_plot as plt2
#-------------------------
# input files
# - field data
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos_TimeSeries.csv'
Litter_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_Litterfall_RawData.csv'

# - lidar canopy profiles
LiDAR_MacHorn_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_LAD_profiles_MacHorn_1m.npz'
LiDAR_rad_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_LAD_profiles_RadiativeTransfer_1m.npz'

# - MODIS data
MODIS_file = '/home/dmilodow/DataStore_DTM/BALI/MODIS/BALI-GEM-MODIS-MOD15A2H-006-results.csv'

# Inventory file, allometry file
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
inventory_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
#-------------------------


# List plots to be assessed and accompanying plot characteristics
plot =        [ 'Belian',   'LF',   'B North','B South',   'E',   'Seraya',   'DC1',    'DC2'  ]
markers =     [    'o',      'o',      'o',      'v',      'v',     'v',       'o',      'v'   ]
edge_colors = [ '#24A23A','#D47E2E','#AC2116','#AC2116','#D47E2E','#24A23A','#24A23A','#24A23A']
fill_colors = [ '#24A23A','#D47E2E','#AC2116','#AC2116','#D47E2E','#24A23A', 'white' , 'white' ]


N = len(plot)
N_sp = 25
subplot_area = 20.*20.
heights = np.arange(0,80.) # for inventory profiles
beta = 0.6 # for inventory profiles
# Some dictionaries to host the data
hemiphot_LAI = {}
LiDAR_MacHorn_LAI = {}
LiDAR_rad_LAI = {}
litter = {}
inventory = {}

# Load in the data
# MODIS
MODIS_LAI = MODIS.load_point_MODIS_LAI_time_series_from_file(MODIS_file, maxQ=True)

# LiDAR
LiDAR1 = np.load(LiDAR_MacHorn_file)
LiDAR2 = np.load(LiDAR_rad_file)

# Inventory
inventory_data = invent.load_crown_survey_data(inventory_file)
a, b, CF, r_sq, p, H, D = invent.retrieve_crown_allometry(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = invent.calculate_allometric_equations_from_survey(inventory_data)

for pp in range(0,N):
    # hemiphotos -> upload subplot level data
    hemiphot_LAI[plot[pp]]={}
    hemiphot_LAI[plot[pp]]['date'],hemiphot_LAI[plot[pp]]['LAI'] = field.get_subplot_LAI_ts(LAI_file,plot[pp])

    # LiDAR data -> retrieve subplot level data
    LiDAR_MacHorn_LAI[plot[pp]] = np.sum(LiDAR1[plot[pp]],axis=1)
    LiDAR_rad_LAI[plot[pp]] = np.sum(LiDAR2[plot[pp]][:,:,-1],axis=1)

    # Inventory file -> retrieve subplot level data
    # set up array to host inventory profiles
    inventory_LAI = np.zeros(N_sp)
    for ss in range(0,N_sp):
        subplot = ss+1
        mask = np.all((inventory_data['plot']==plot[pp],inventory_data['subplot']==subplot),axis=0)
        Ht,Area,Depth = invent.calculate_crown_dimensions(inventory_data['DBH_field'][mask],inventory_data['Height_field'][mask],inventory_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        inventory_LAD, CanopyV = invent.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)
        inventory_LAI[ss] = np.sum(inventory_LAD)
    inventory['volume'] = inventory_LAI.copy()
    

    # Litter file -> upload subplot level data
    litter[plot[pp]]={}
    litter[plot[pp]]['date'],litter[plot[pp]]['previous date'],litter[plot[pp]]['flux'] = field.get_subplot_litterfall_ts(Litter_file,plot[pp])


#------------------------------------------------------------------------------------------------------------
# Now plot up the data.  First - time series of LAI observations
plt.figure(1, facecolor='White',figsize=[9,6])

ax1 = plt.subplot2grid((3,1),(0,0))
ax2 = plt.subplot2grid((3,1),(1,0),sharex=ax1,sharey=ax1)
ax3 = plt.subplot2grid((3,1),(2,0),sharex=ax1)
for pp in range(0,N):
    MODIS_plot = plot[pp]
    if MODIS_plot == 'B North':
        MODIS_plot = 'BNorth'
    if MODIS_plot == 'B South':
        MODIS_plot = 'BSouth'

    ax1.plot(MODIS_LAI[MODIS_plot]['date'].astype(datetime),MODIS_LAI[MODIS_plot]['LAI'],ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])

    ax2.plot(hemiphot_LAI[plot[pp]]['date'].astype(datetime),np.mean(hemiphot_LAI[plot[pp]]['LAI'],axis=0),ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


    ax3.plot(litter[plot[pp]]['date'].astype(datetime),np.mean(litter[plot[pp]]['flux'],axis=0),ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


# Configure plots
ax1.annotate('a - MODIS LAI', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
ax1.set_ylabel('LAI',fontsize=axis_size)
ax2.annotate('b - Hemisfer LAI', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
ax2.set_ylabel('LAI',fontsize=axis_size)
ax3.annotate('a - Litter traps', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
ax3.set_ylabel('litter flux / Mg ha$^{-1}$ yr$^{-1}$',fontsize=axis_size)

ax3.legend(loc='center right',fontsize = rcParams['font.size'])
ax1.set_xlim(xmin=np.datetime64('2011-06-01','D').astype(datetime))
plt.tight_layout()
plt.show()

