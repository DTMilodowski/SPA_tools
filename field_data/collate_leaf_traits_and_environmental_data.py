import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import sys

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/')
#import statistics_tools as stats2
import load_field_data as field
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/')
import LiDAR_tools as lidar
import inventory_based_LAD_profiles as inventory
import LiDAR_MacHorn_LAD_profiles as LAD
import auxilliary_functions as aux
import canopy_microclimate as clim
from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2



chem_file = '../../../BALI_traits_data/CombinedPlots/BALI_traits_CN_29032017.csv'
spp_file = '../../../BALI_traits_data/CombinedPlots/BALI_species_28022017.csv'
photo_file = '../../../BALI_traits_data/CombinedPlots/Photosynthesis_Combined_14122016_hyphens.csv'
branch_file = '../../../BALI_traits_data/CombinedPlots/ParameterTrees_Combined_14122016.csv'
leaf_file = '../../../BALI_traits_data/CombinedPlots/LeafArea_Combined_14122016.csv'
census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
light_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/src/output/BALI_subplot_lighttransmittance.npz'
subplot_coordinate_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/BALI_subplot_coordinates_corrected.csv' #check
las_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/Carbon_plot_point_cloud_buffer.las'

# get LiDAR data and sort into plots
lidar_pts = lidar.load_lidar_data(las_file)
plot_pts = {}

# get lower left hand corner of plots
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
plot_origin = {}
buffer_width = 5.
for i in range(0,len(subplot_polygons.keys())):
    plot = subplot_polygons.keys()[i]
    print plot
    plot_origin[plot] = subplot_polygons[plot][0,0,:]

    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[plot].shape[0]*subplot_polygons[plot].shape[1]
    coord_pairs = subplot_polygons[plot].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box_with_buffer(coord_pairs,buffer_width)
    plot_pts[plot] = lidar.filter_lidar_data_by_polygon(lidar_pts,bbox_polygon)

# clear up some memory
lidar_pts = None


# Load in the traits
branch, spp, genus, N, Narea, C, CNratio, SLA, LMA, LeafArea, LeafThickness, LeafHeight, VPD, Rd, Vcmax, Jmax, ShadeTag, ftype = field.collate_branch_level_traits(chem_file,photo_file,branch_file,leaf_file,spp_file)

# Filter traits data
LeafThickness[LeafThickness>0.60]=np.nan
Vcmax[Vcmax>150]=np.nan
#LMA[LMA>0.025]=np.nan
N[N>4]=np.nan
#Narea*=10.**3
Narea[Narea>4.0]=np.nan
LMA_C = LMA*C/100. # convert to g(C) m^-2
LMA_C[LMA_C>100.]=np.nan

# Loop through branches pulling out plot codes
N_branches = branch.size
plot_list = []
tree_list = []
for i in range(0,N_branches):
    plot_temp = branch[i].split('-')[0]
    if plot_temp =='BEL':
        plot_list.append('Belian')
    elif plot_temp =='SER':
        plot_list.append('Seraya')
    elif plot_temp =='DAS1':
        plot_list.append('DC1')
    elif plot_temp =='DAF2':
        plot_list.append('DC2')
    elif plot_temp =='ESA':
        plot_list.append('E')
    elif plot_temp =='BNT':
        plot_list.append('B North')
    elif plot_temp =='BSO':
        plot_list.append('B South')
    elif plot_temp =='SLF':
        plot_list.append('LF')
    else:
        print plot_temp
    if branch[i].split('-')[1][1]=='B': # find out why this tree labelled as such
        tree_list.append(float(branch[i].split('-')[1][2:]))
    else:
        tree_list.append(float(branch[i].split('-')[1][1:]))

# Now need to find subplot in which this tree is located
plot = np.asarray(plot_list)
tree = np.asarray(tree_list)

census_plot, census_subplot, census_dates, tree_tag, alt_tag, DPOM, HPOM, TreeHeight, C_stem, C_coarse_root, RAINFOR, Alive_flag, census_spp, SubplotCoords, WoodDensity = field.read_ICP_census_data(census_file)
field_data = inventory.load_crown_survey_data(field_file)


subplot = np.zeros(N_branches)*np.nan
light_availability = np.zeros(N_branches)*np.nan
tree_centric = np.zeros(N_branches) # a flag to state whether or not we have tree-centric LAD to estimate light environment

plots = np.unique(plot)
radius = 5
layer_thickness = 1
max_height = 80
k=0.5
for i in range(0,N_branches):
    tree_index = np.all((census_plot==plot[i], np.any((tree_tag == tree[i], alt_tag == tree[i]),axis=0)),axis=0)
    tree_index_alt = np.all((field_data['plot']==plot[i], field_data['tag'] == tree[i]),axis=0)

    # first case scenario - we have located tree in canopy map
    if tree_index_alt.sum() == 1:
        subplot[i] = field_data['subplot'][tree_index_alt]
        tree_xy = np.asarray([field_data['Xfield'][tree_index_alt]+ plot_origin[plot[i]][0],field_data['Yfield'][tree_index_alt]+ plot_origin[plot[i]][1]])
        tree_pts = filter_lidar_data_by_neighbourhood(plot_pts[plot[i]],tree_xy,radius)
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(tree_pts, max_height,layer_thickness)
        LAD_tree = LAD.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)
        I_tree=clim.estimate_canopy_light_transmittance(LAD_tree,heights,k)
        ht_index = int(LeafHeight[i])
        light_availability[i]=I_tree[ht_index]
        tree_centric[i] = 1

    # otherwise, need to use census data, which only gives subplot level data
    elif tree_index.sum() == 1:
        #print "found tree :-)"
        subplot[i] = census_subplot[tree_index]
        subplot_poly = subplot_polygons[Plot_name][subplot_labels[plot[i]] == subplot[i],:,:]
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_pts[plot[i]],subplot_poly)
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(tree_pts, max_height,layer_thickness)
        LAD_subplot = LAD.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)
        I_subplot=clim.estimate_canopy_light_transmittance(LAD_subplot,heights,k)
        ht_index = int(LeafHeight[i])
        light_availability[i]=I_subplot[ht_index]

    elif np.all((tree_index.sum()==0,tree_index_alt.sum()==0)):
        print "cannot find tree from traits record in the census data"
        print "\t tree: ", tree[i], "; plot: ", plot[i]
        subplot[i] = np.nan
        light_availability[i]=np.nan
    else:
        print "we have a problem - more than one tree in this plot has the same tag "
        subplot[i] = np.nan         
        light_availability[i]=np.nan


# write traits data to a csv file for ingestion into R
#plot,subplot,forest_type,branch,shade_tag,spp,genus,leaf_height,light_availability,leaf_thickness,leaf_area,LMA,C%,Carea,N%,Narea,Vcmax,Rd
out = open('BALI_leaf_traits.csv','w')
out.write('plot,subplot,forest_type,branch,shade_tag,spp,genus,leaf_height,light_availability,leaf_thickness,leaf_area,LMA,C%,Carea,N%,Narea,Vcmax,Rd\n')
for i in range(0,N_branches):
    out.write(plot[i] + ',' + str(subplot[i])  + ', ' + ftype[i] + ',' + branch[i]  + ',' + str(ShadeTag[i])  + ',' + spp[i] + ',' + spp[i].split(' ')[0] + ',' + str(LeafHeight[i]) + ',' + str(light_availability[i]) + + ',' + str(tree_centric[i]) + ',' +str( LeafThickness[i]) + ',' +str( LeafArea[i]) + ',' + str(LMA[i]) + ',' + str(C[i]) + ',' + str(LMA_C[i]) + ',' + str(N[i]) + ',' + str(Narea[i]) + ',' + str(Vcmax[i]) + ',' + str(Rd[i]) + '\n')

out.close()

