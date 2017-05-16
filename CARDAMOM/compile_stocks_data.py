# This function essentially compiles together the carbon stocks data and plots up time
# series - both at subplot level and a compilation plot at the 1ha level.
#----------------------------------------------------------------------------------------
import numpy as np
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/')
import load_field_data as field
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
census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'

# List plots to be assessed and accompanying plot characteristics
plot =        [ 'Belian',   'Seraya',  'LF',     'E'   ,'B North','B South']
markers =     [    'o',       'v'  ,   'o',      'v'   ,   'o',      'v'   ]
edge_colors = [ '#24A23A','#24A23A','#D47E2E','#D47E2E','#AC2116','#AC2116']
fill_colors = [ 'white','#24A23A','white','#D47E2E','white','#AC2116']


N = len(plot)
N_sp = 25

# Load in the data
# Census data
census = field.collate_plot_level_census_data(census_file)


# plot census data - use a multipanel plot showing individual subplot evolution (1 per plot) and compiled 1ha plots
plt.figure(1, facecolor='White',figsize=[9,9])
ax1 = plt.subplot2grid((3,3),(0,0))
ax2 = plt.subplot2grid((3,3),(0,1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3,3),(0,2),sharey=ax1,sharex=ax1)
ax4 = plt.subplot2grid((3,3),(1,0),sharey=ax1,sharex=ax1)
ax5 = plt.subplot2grid((3,3),(1,1),sharey=ax1,sharex=ax1)
ax6 = plt.subplot2grid((3,3),(1,2),sharey=ax1,sharex=ax1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=4, sharex=ax1)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]
ann = ['a','b','c','d','e','f']

# first up lets deal with the subplots
for pp in range(0,N):
    axes[pp].annotate(ann[pp]+' - '+ plot[pp], xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
    for ss in range(0,N_sp):
        Cwood = census[plot[pp]]['C_wood'][ss,:]*25./1000.
        date = census[plot[pp]]['CensusDate'][ss,:]
        Cwood[date<np.datetime64('2000-01-01','D')] = np.nan
        axes[pp].plot(date.astype(datetime),Cwood,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp])

    # now lets plot the averages
    Cwood = np.sum(census[plot[pp]]['C_wood'],axis=0)/1000.
    date = np.max(census[plot[pp]]['CensusDate'],axis=0)
    Cwood[date<np.datetime64('2000-01-01','D')] = np.nan
    ax7.plot(date.astype(datetime),Cwood,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


# configure plot
ax7.annotate('g - all SAFE 1 ha plots', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() +ax3.get_xticklabels()
yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() +ax5.get_yticklabels() +ax6.get_yticklabels()

ax1.set_ylabel('Cwood / Mg ha$^{-1}$', fontsize = axis_size)
ax4.set_ylabel('Cwood / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylabel('Cwood / Mg ha$^{-1}$', fontsize = axis_size)
ax6.set_ylim(ymax=1000)
ax7.set_ylim(ymax=500)
ax7.set_xlim(xmax=np.datetime64('2016-01-01').astype(datetime))
ax7.legend(loc='upper right',fontsize = rcParams['font.size'])
plt.setp(xticklabels,visible=False)
plt.setp(yticklabels,visible=False)

plt.subplots_adjust(hspace=0.001,wspace=0.001)
plt.tight_layout()
plt.savefig('SAFE_Cwood_stocks.png')
plt.show()

# Now plot growth
plt.figure(2, facecolor='White',figsize=[9,9])
ax1 = plt.subplot2grid((3,3),(0,0))
ax2 = plt.subplot2grid((3,3),(0,1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3,3),(0,2),sharey=ax1,sharex=ax1)
ax4 = plt.subplot2grid((3,3),(1,0),sharey=ax1,sharex=ax1)
ax5 = plt.subplot2grid((3,3),(1,1),sharey=ax1,sharex=ax1)
ax6 = plt.subplot2grid((3,3),(1,2),sharey=ax1,sharex=ax1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=4, sharex=ax1)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]
ann = ['a','b','c','d','e','f']

# first up lets deal with the subplots
for pp in range(0,N):
    axes[pp].annotate(ann[pp]+' - '+ plot[pp], xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
    growth = census[plot[pp]]['Growth']/1000.
    growth[np.isnan(growth)]=0.
    growth= np.cumsum(growth,axis=1)
    for ss in range(0,N_sp):
        sp_growth= growth[ss,:]*25.
        date = census[plot[pp]]['CensusDate'][ss,:]
        sp_growth[date<np.datetime64('2000-01-01','D')] = np.nan
        axes[pp].plot(date.astype(datetime),sp_growth,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp])

    # now lets plot the cumulative sum
    plot_growth = np.sum(growth,axis=0)
    date = np.max(census[plot[pp]]['CensusDate'],axis=0)
    plot_growth[date<np.datetime64('2000-01-01','D')] = np.nan
    ax7.plot(date.astype(datetime),plot_growth,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])

# configure plot
ax7.annotate('g - all SAFE 1 ha plots', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() +ax3.get_xticklabels()
yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() +ax5.get_yticklabels() +ax6.get_yticklabels()

ax1.set_ylabel('growth / Mg ha$^{-1}$', fontsize = axis_size)
ax4.set_ylabel('growth / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylabel('growth / Mg ha$^{-1}$', fontsize = axis_size)
#ax7.set_ylim(ymax=400000)
ax7.set_xlim(xmax=np.datetime64('2016-01-01').astype(datetime))
ax7.legend(loc='center left',fontsize = rcParams['font.size'])
plt.setp(xticklabels,visible=False)
plt.setp(yticklabels,visible=False)

plt.subplots_adjust(hspace=0.001,wspace=0.001)
plt.tight_layout()
plt.savefig('SAFE_Cwood_growth.png')
plt.show()




# Now plot recruitment
plt.figure(3, facecolor='White',figsize=[9,9])
ax1 = plt.subplot2grid((3,3),(0,0))
ax2 = plt.subplot2grid((3,3),(0,1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3,3),(0,2),sharey=ax1,sharex=ax1)
ax4 = plt.subplot2grid((3,3),(1,0),sharey=ax1,sharex=ax1)
ax5 = plt.subplot2grid((3,3),(1,1),sharey=ax1,sharex=ax1)
ax6 = plt.subplot2grid((3,3),(1,2),sharey=ax1,sharex=ax1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=4, sharex=ax1)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]
ann = ['a','b','c','d','e','f']

# first up lets deal with the subplots
for pp in range(0,N):
    axes[pp].annotate(ann[pp]+' - '+ plot[pp], xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
    recruit = census[plot[pp]]['Recruitment']/1000.
    recruit[np.isnan(recruit)]=0.
    recruit= np.cumsum(recruit,axis=1)
    for ss in range(0,N_sp):
        sp_recruit= recruit[ss,:]*25.
        date = census[plot[pp]]['CensusDate'][ss,:]
        sp_recruit[date<np.datetime64('2000-01-01','D')] = np.nan
        axes[pp].plot(date.astype(datetime),sp_recruit,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp])

    # now lets plot the cumulative sum
    plot_recruit = np.sum(recruit,axis=0)
    date = np.max(census[plot[pp]]['CensusDate'],axis=0)
    plot_recruit[date<np.datetime64('2000-01-01','D')] = np.nan
    ax7.plot(date.astype(datetime),plot_recruit,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


# configure plot
ax7.annotate('g - all SAFE 1 ha plots', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() +ax3.get_xticklabels()
yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() +ax5.get_yticklabels() +ax6.get_yticklabels()

ax1.set_ylabel('recruitment / Mg ha$^{-1}$', fontsize = axis_size)
ax4.set_ylabel('recruitment / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylabel('recruitment / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylim(ymax=10)
ax7.set_xlim(xmax=np.datetime64('2016-01-01').astype(datetime))
ax7.legend(loc='center left',fontsize = rcParams['font.size'])
plt.setp(xticklabels,visible=False)
plt.setp(yticklabels,visible=False)

plt.subplots_adjust(hspace=0.001,wspace=0.001)
plt.tight_layout()
plt.savefig('SAFE_Cwood_recruitment.png')
plt.show()


# Now plot mortality
plt.figure(4, facecolor='White',figsize=[9,9])
ax1 = plt.subplot2grid((3,3),(0,0))
ax2 = plt.subplot2grid((3,3),(0,1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3,3),(0,2),sharey=ax1,sharex=ax1)
ax4 = plt.subplot2grid((3,3),(1,0),sharey=ax1,sharex=ax1)
ax5 = plt.subplot2grid((3,3),(1,1),sharey=ax1,sharex=ax1)
ax6 = plt.subplot2grid((3,3),(1,2),sharey=ax1,sharex=ax1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=4, sharex=ax1)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]
ann = ['a','b','c','d','e','f']

# first up lets deal with the subplots
for pp in range(0,N):
    axes[pp].annotate(ann[pp]+' - '+ plot[pp], xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
    mort = census[plot[pp]]['Mortality']/1000.
    mort[np.isnan(mort)]=0.
    mort= np.cumsum(mort,axis=1)
    for ss in range(0,N_sp):
        sp_mort= mort[ss,:]*25.
        date = census[plot[pp]]['CensusDate'][ss,:]
        sp_mort[date<np.datetime64('2000-01-01','D')] = np.nan
        axes[pp].plot(date.astype(datetime),sp_mort,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp])

    # now lets plot the cumulative sum
    plot_mort = np.sum(mort,axis=0)
    date = np.max(census[plot[pp]]['CensusDate'],axis=0)
    plot_mort[date<np.datetime64('2000-01-01','D')] = np.nan
    ax7.plot(date.astype(datetime),plot_mort,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


# configure plot
ax7.annotate('g - all SAFE 1 ha plots', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() +ax3.get_xticklabels()
yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() +ax5.get_yticklabels() +ax6.get_yticklabels()

ax1.set_ylabel('mortality / Mg ha$^{-1}$', fontsize = axis_size)
ax4.set_ylabel('mortality / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylabel('mortality / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_xlim(xmax=np.datetime64('2016-01-01').astype(datetime))
ax7.legend(loc='center left',fontsize = rcParams['font.size'])
plt.setp(xticklabels,visible=False)
plt.setp(yticklabels,visible=False)

plt.subplots_adjust(hspace=0.001,wspace=0.001)
plt.tight_layout()
plt.savefig('SAFE_Cwood_mortality.png')
plt.show()


# Now plot NPP
plt.figure(5, facecolor='White',figsize=[9,9])
ax1 = plt.subplot2grid((3,3),(0,0))
ax2 = plt.subplot2grid((3,3),(0,1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3,3),(0,2),sharey=ax1,sharex=ax1)
ax4 = plt.subplot2grid((3,3),(1,0),sharey=ax1,sharex=ax1)
ax5 = plt.subplot2grid((3,3),(1,1),sharey=ax1,sharex=ax1)
ax6 = plt.subplot2grid((3,3),(1,2),sharey=ax1,sharex=ax1)
ax7 = plt.subplot2grid((3,3),(2,0),colspan=4, sharex=ax1)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]
ann = ['a','b','c','d','e','f']

# first up lets deal with the subplots
for pp in range(0,N):
    axes[pp].annotate(ann[pp]+' - '+ plot[pp], xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 
    npp = (census[plot[pp]]['Growth'] + census[plot[pp]]['Recruitment'] -  census[plot[pp]]['Mortality'])/1000.
    npp[np.isnan(npp)]=0.
    npp= np.cumsum(npp,axis=1)
    for ss in range(0,N_sp):
        sp_npp= npp[ss,:]*25.
        date = census[plot[pp]]['CensusDate'][ss,:]
        sp_npp[date<np.datetime64('2000-01-01','D')] = np.nan
        axes[pp].plot(date.astype(datetime),sp_npp,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp])

    # now lets plot the cumulative sum
    plot_npp = np.sum(npp,axis=0)
    date = np.max(census[plot[pp]]['CensusDate'],axis=0)
    plot_npp[date<np.datetime64('2000-01-01','D')] = np.nan
    ax7.plot(date.astype(datetime),plot_npp,ls=':',c=edge_colors[pp],marker=markers[pp],mec=edge_colors[pp],mfc=fill_colors[pp],label=plot[pp])


# configure plot
ax7.annotate('g - all SAFE 1 ha plots', xy=(0.02,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2,backgroundcolor='white') 

xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() +ax3.get_xticklabels()
yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels() +ax5.get_yticklabels() +ax6.get_yticklabels()

ax1.set_ylabel('cumulative NPP / Mg ha$^{-1}$', fontsize = axis_size)
ax4.set_ylabel('cumulative NPP / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_ylabel('cumulative NPP / Mg ha$^{-1}$', fontsize = axis_size)
ax7.set_xlim(xmax=np.datetime64('2016-01-01').astype(datetime))
ax7.legend(loc='center left',fontsize = rcParams['font.size'])
plt.setp(xticklabels,visible=False)
plt.setp(yticklabels,visible=False)

plt.subplots_adjust(hspace=0.001,wspace=0.001)
plt.tight_layout()
plt.savefig('SAFE_Cwood_npp.png')
plt.show()

