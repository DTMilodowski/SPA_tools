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
    if plot == 'Seraya':
        plot_origin[plot] = subplot_polygons[plot][subplot_labels['Seraya']==6,0,:][0]
    elif plot == 'DC1':
        plot_origin[plot] = subplot_polygons[plot][subplot_labels['DC1']==1,3,:][0]
    elif plot == 'DC2':
        plot_origin[plot] = subplot_polygons[plot][subplot_labels['DC2']==1,3,:][0]
    else:
        plot_origin[plot] = subplot_polygons[plot][0,0,:]

    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[plot].shape[0]*subplot_polygons[plot].shape[1]
    coord_pairs = subplot_polygons[plot].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box_with_buffer(coord_pairs,buffer_width)
    plot_pts[plot] = lidar.filter_lidar_data_by_polygon(lidar_pts,bbox_polygon)

# clear up some memory
lidar_pts = None


# Load in the traits
cleaning_scheme = 1 # 0 = mine, 1 = Sabine's, 2 = conservative
branch, spp, genus, N, Narea, C, CNratio, SLA, LMA, LeafArea, LeafThickness, LeafHeight, VPD, Rd, Vcmax, Jmax, ShadeTag, gs, ftype = field.collate_branch_level_traits(chem_file,photo_file,branch_file,leaf_file,spp_file,cleaning_scheme)

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
radius = 10.
layer_thickness = 1
max_height = 80
k=0.5
for i in range(0,N_branches):
    print i+1, '/', N_branches, plot[i]
    tree_index = np.all((census_plot==plot[i], np.any((tree_tag == tree[i], alt_tag == tree[i]),axis=0)),axis=0)
    tree_index_alt = np.all((field_data['plot']==plot[i], field_data['tag'] == tree[i]),axis=0)

    # first case scenario - we have located tree in canopy map
    if tree_index_alt.sum() == 1:
        subplot[i] = field_data['subplot'][tree_index_alt]
        tree_xy = np.asarray([field_data['Xfield'][tree_index_alt]+ plot_origin[plot[i]][0],field_data['Yfield'][tree_index_alt]+ plot_origin[plot[i]][1]])
        tree_pts = lidar.filter_lidar_data_by_neighbourhood(plot_pts[plot[i]],tree_xy,radius)
        if tree_pts.size > 0:
            heights,first_return_profile,n_ground_returns = LAD.bin_returns(tree_pts, max_height,layer_thickness)
            LAD_tree = LAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)
            I_tree=clim.estimate_canopy_light_transmittance(LAD_tree,heights,k)
            """
            # subplot based profile
            subplot_poly = subplot_polygons[plot[i]][subplot_labels[plot[i]] == subplot[i],:,:][0]
            sp_pts = lidar.filter_lidar_data_by_polygon(plot_pts[plot[i]],subplot_poly)
            heights_sp,first_return_profile_sp,n_ground_returns_sp = LAD.bin_returns(sp_pts, max_height,layer_thickness)
            LAD_sp = LAD.estimate_LAD_MacArthurHorn(first_return_profile_sp, n_ground_returns_sp, layer_thickness, 1.)
            I_sp=clim.estimate_canopy_light_transmittance(LAD_sp,heights_sp,k)
            """
            if np.isfinite(LeafHeight[i]):
                ht_index = int(LeafHeight[i])
                light_availability[i]=I_tree[ht_index]
                tree_centric[i] = 1
            else:
                light_availability[i]=np.nan
            #plt.show()
        else:
            light_availability[i]=np.nan
            

    # otherwise, need to use census data, which only gives subplot level data
    elif tree_index.sum() == 1:
        #print "found tree :-)"
        subplot[i] = census_subplot[tree_index]
        subplot_poly = subplot_polygons[plot[i]][subplot_labels[plot[i]] == subplot[i],:,:]
        tree_pts = lidar.filter_lidar_data_by_polygon(plot_pts[plot[i]],subplot_poly)
        if tree_pts.size > 0:
            heights,first_return_profile,n_ground_returns = LAD1.bin_returns(tree_pts, max_height,layer_thickness)
            LAD_subplot = LAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)
            I_subplot=clim.estimate_canopy_light_transmittance(LAD_subplot,heights,k)
            if np.isfinite(LeafHeight[i]):
                ht_index = int(LeafHeight[i])
                light_availability[i]=I_subplot[ht_index]
            else:
                light_availability[i]=I_tree[ht_index]=np.nan
        else:
            light_availability[i]=I_tree[ht_index]=np.nan

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
out.write('plot,subplot,forest_type,branch,shade_tag,spp,genus,leaf_height,light_availability,tree_centric,leaf_thickness,leaf_area,LMA,C%,Carea,N%,Narea,Vcmax,Rd\n')
for i in range(0,N_branches):
    out.write(plot[i] + ',' + str(subplot[i])  + ', ' + ftype[i] + ',' + branch[i]  + ',' + str(ShadeTag[i])  + ',' + spp[i] + ',' + spp[i].split(' ')[0] + ',' + str(LeafHeight[i]) + ',' + str(light_availability[i]) + ',' + str(tree_centric[i]) + ',' + str( LeafThickness[i]) + ',' + str( LeafArea[i]) + ',' + str(LMA[i]) + ',' + str(C[i]) + ',' + str(LMA_C[i]) + ',' + str(N[i]) + ',' + str(Narea[i]) + ',' + str(Vcmax[i]) + ',' + str(Rd[i]) + '\n')

out.close()



"""
# let's do some checks on positioning
for i in range(0,len(plots)):
    
    plt.figure(1, facecolor='White',figsize=[6,6])
    ax1 = plt.subplot2grid((1,1),(0,0))
    print plots[i]
    for s in range(0,25):
        plt.plot(subplot_polygons[plots[i]][s,:,0],subplot_polygons[plots[i]][s,:,1],'-')
    plt.plot(plot_origin[plots[i]][0],plot_origin[plots[i]][1],'o')

    mask = field_data['plot']==plots[i]
    tree_xy = np.asarray([field_data['Xfield'][mask]+ plot_origin[plots[i]][0],field_data['Yfield'][mask]+ plot_origin[plots[i]][1]])
    #print tree_xy
    plt.plot(tree_xy[0,:],tree_xy[1,:],'o')
    plt.show()

"""

"""

# Figure 4 - vertical distribution of leaf traits in terms of light availability rather than leaf height per se:
# Vcmax, Rd, N,
# LMA, Leaf Thickness, C:N ratio
plt.figure(4, facecolor='White',figsize=[9,6])

ax1 = plt.subplot2grid((2,3),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Vcmax[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Vcmax[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Vcmax) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Vcmax[mask],light_availability[mask])

ax1.plot(Vcmax[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax1.plot(Vcmax[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax1.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1.set_ylabel('Light availability')
ax1.set_xlabel('$V_{cmax}$')

ax2 = plt.subplot2grid((2,3),(0,1),sharey=ax1)
ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(Rd[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax2.plot(Rd[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Rd[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Rd[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Rd[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Rd[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Rd) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Rd[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '


pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_ylabel('Light availability')
ax2.set_xlabel('$R_d$')

ax3 = plt.subplot2grid((2,3),(0,2),sharey=ax1)
ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
#ax3.plot(N[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
#ax3.plot(N[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')
ax3.plot(Narea[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax3.plot(Narea[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

#mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
#m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
#mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
#m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Narea) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Narea[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_ylabel('Light availability')
#ax3.set_xlabel('%N')
ax3.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')


ax4 = plt.subplot2grid((2,3),(1,0),sharey=ax1)
ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax4.plot(LMA_C[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax4.plot(LMA_C[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')
print "LMA(C)"
mask = ~np.isnan(LMA_C[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
print "OG: ", np.mean(LMA_C[ftype=='OG'][mask]), " +/- ", np.std(LMA_C[ftype=='OG'][mask])
m1, c1, r1, p1, err1 = stats.linregress(LMA[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(LMA_C[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
print "SL: ", np.mean(LMA_C[ftype=='SL'][mask]), " +/- ", np.std(LMA_C[ftype=='SL'][mask])
m2, c2, r2, p2, err2 = stats.linregress(LMA[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(LMA_C) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(LMA_C[mask],light_availability[mask])

mask = ~np.isnan(LMA_C) & ~np.isnan(light_availability)
print "all: ", np.mean(LMA_C[mask]), " +/- ", np.std(LMA_C[mask])

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax4.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax4.set_ylabel('Light availability')
ax4.set_xlabel('LMA(C) / g m$^{-2}$')

ax5 = plt.subplot2grid((2,3),(1,1),sharey=ax1)
ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax5.plot(LeafThickness[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax5.plot(LeafThickness[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

mask = ~np.isnan(LeafThickness[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(LeafThickness[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(LeafThickness[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(LeafThickness[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(LeafThickness) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(LeafThickness[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax5.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax5.set_ylabel('Light availability')
ax5.set_xlabel('Leaf Thickness / mm')

ax6 = plt.subplot2grid((2,3),(1,2),sharey=ax1)
ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax6.plot(CNratio[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax6.plot(CNratio[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

mask = ~np.isnan(CNratio[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(CNratio[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(CNratio[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(CNratio[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(CNratio) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(CNratio[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax6.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax6.set_ylabel('Light availability')
ax6.set_xlabel('C:N ratio')

plt.tight_layout()
plt.savefig('vertical_trait_distributions_with_respect_to_light.png')


# plot light vs height
plt.figure(5, facecolor='White',figsize=[5,5])

ax1 = plt.subplot2grid((1,1),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(light_availability[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(light_availability[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(light_availability[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(light_availability[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])

ax1.plot(light_availability[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax1.plot(light_availability[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str
ax1.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1.set_ylabel('Height / m')
ax1.set_xlabel('light_availability')
plt.savefig('height_vs_light.png')
plt.show()


"""


# Figure 6 - vertical distribution of chemo-photosynthetic leaf traits in terms of light availability and leaf height, for comparison:
# Vcmax, Rd, N,
plt.figure(6, facecolor='White',figsize=[6,9])

# axis 2 - height vs. Vcmax
ax1 = plt.subplot2grid((3,2),(1,0))
ax1.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Vcmax[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Vcmax[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
mask = ~np.isnan(Vcmax) & ~np.isnan(LeafHeight)
m, c, r, p, err = stats.linregress(Vcmax[mask],LeafHeight[mask])

ax1.plot(Vcmax[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax1.plot(Vcmax[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax1.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1.set_ylabel('Leaf Height / m')
ax1.set_xlabel('$V_{cmax}$')


# axis 2 - light availability vs. Vcmax
ax1b = plt.subplot2grid((3,2),(1,1))
ax1b.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Vcmax[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Vcmax[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Vcmax) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Vcmax[mask],light_availability[mask])

ax1b.plot(Vcmax[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax1b.plot(Vcmax[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '



stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax1b.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1b.set_ylabel('Light availability')
ax1b.set_xlabel('$V_{cmax}$')


# Axis 2 Rd vs height
ax2 = plt.subplot2grid((3,2),(2,0),sharey=ax1)
ax2.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(Rd[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax2.plot(Rd[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Rd[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Rd[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Rd[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Rd[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
mask = ~np.isnan(Rd) & ~np.isnan(LeafHeight)
m, c, r, p, err = stats.linregress(Rd[mask],LeafHeight[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '


pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_ylabel('Leaf Height / m')
ax2.set_xlabel('$R_d$')

# Axis 2b Rd vs light availability
ax2b = plt.subplot2grid((3,2),(2,1),sharey=ax1b)
ax2b.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2b.plot(Rd[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax2b.plot(Rd[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Rd[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Rd[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Rd[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Rd[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Rd) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Rd[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '


pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax2b.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2b.set_ylabel('Light availability')
ax2b.set_xlabel('$R_d$')


# Axis 3 height vs Narea
ax3 = plt.subplot2grid((3,2),(0,0),sharey=ax1)
ax3.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax3.plot(Narea[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax3.plot(Narea[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
mask = ~np.isnan(Narea) & ~np.isnan(LeafHeight)
m, c, r, p, err = stats.linregress(Narea[mask],LeafHeight[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '
                                   
stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_ylabel('Leaf Height / m')
ax3.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')


# N-area vs light
ax3b = plt.subplot2grid((3,2),(0,1),sharey=ax1b)
ax3b.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax3b.plot(Narea[ftype=='OG'],light_availability[ftype=='OG'],'.',color='blue')
ax3b.plot(Narea[ftype=='SL'],light_availability[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(light_availability[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],light_availability[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(light_availability[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],light_availability[ftype=='SL'][mask])
mask = ~np.isnan(Narea) & ~np.isnan(light_availability)
m, c, r, p, err = stats.linregress(Narea[mask],light_availability[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '
                                     
stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax3b.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3b.set_ylabel('Light availability')
ax3b.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')

plt.tight_layout()
plt.savefig('vertical_chemophoto_trait_distributions.png')



# Figure 7 - photosynthesis vs. N:
# Vcmax, Rd, N,
plt.figure(7, facecolor='White',figsize=[6,6])

# axis 2 - Narea vs. Vcmax
ax1 = plt.subplot2grid((2,2),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(Narea[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],Vcmax[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],Vcmax[ftype=='SL'][mask])
mask = ~np.isnan(Vcmax) & ~np.isnan(Narea)
m, c, r, p, err = stats.linregress(Narea[mask],Vcmax[mask])

ax1.plot(Narea[ftype=='OG'],Vcmax[ftype=='OG'],'.',color='blue')
ax1.plot(Narea[ftype=='SL'],Vcmax[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax1.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')
ax1.set_ylabel('$V_{cmax}$')

# Axis 2 Rd vs Narea
ax2 = plt.subplot2grid((2,2),(0,1),sharex=ax1)
ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(Narea[ftype=='OG'],Rd[ftype=='OG'],'.',color='blue')
ax2.plot(Narea[ftype=='SL'],Rd[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Rd[ftype=='OG']) & ~np.isnan(Narea[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],Rd[ftype=='OG'][mask])
mask = ~np.isnan(Rd[ftype=='SL']) & ~np.isnan(Narea[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],Rd[ftype=='SL'][mask])
mask = ~np.isnan(Rd) & ~np.isnan(Narea)
m, c, r, p, err = stats.linregress(Narea[mask],Rd[mask])
p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')
ax2.set_ylabel('$R_d$')

# axis 3 - Stomatal conductance vs. Vcmax
ax3 = plt.subplot2grid((2,2),(1,0),sharey=ax1)
ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(gs[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(gs[ftype=='OG'][mask],Vcmax[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(gs[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(gs[ftype=='SL'][mask],Vcmax[ftype=='SL'][mask])
mask = ~np.isnan(Vcmax) & ~np.isnan(gs)
m, c, r, p, err = stats.linregress(gs[mask],Vcmax[mask])

ax3.plot(gs[ftype=='OG'],Vcmax[ftype=='OG'],'.',color='blue')
ax3.plot(gs[ftype=='SL'],Vcmax[ftype=='SL'],'.',color='red')

p1str=''
if p1<0.001:
    p1str='***'
elif p1<0.01:
    p1str='** '
elif p1<0.05:
    p1str='*  '
elif p1<0.1:
    p1str='$^.$  '
else:
    p1str='   '

p2str=''
if p2<0.001:
    p2str='***'
elif p2<0.01:
    p2str='** '
elif p2<0.05:
    p2str='*  '
elif p2<0.1:
    p2str='$^.$  '
else:
    p2str='   '

pstr=''
if p<0.001:
    pstr='***'
elif p<0.01:
    pstr='** '
elif p<0.05:
    pstr='*  '
elif p<0.1:
    pstr='$^.$  '
else:
    pstr='   '

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str + '\nR$^2$=' + '%.3f' % r**2 + pstr
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_xlabel('stomatal conductance')
ax3.set_ylabel('$V_{cmax}$')



plt.tight_layout()
plt.savefig('photo_N.png')
plt.show()