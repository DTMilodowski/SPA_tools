import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import sys

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/')
#import statistics_tools as stats2
import load_field_data as field

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

# Load in the traits
branch, spp, genus, N, Narea, C, CNratio, SLA, LMA, LeafArea, LeafThickness, LeafHeight, VPD, Rd, Vcmax, Jmax, ShadeTag, ftype = field.collate_branch_level_traits(chem_file,photo_file,branch_file,leaf_file,spp_file)

LeafThickness[LeafThickness>0.60]=np.nan
Vcmax[Vcmax>150]=np.nan
#LMA[LMA>0.025]=np.nan
N[N>4]=np.nan
#Narea*=10.**3
Narea[Narea>4.0]=np.nan
LMA_C = LMA*C/100. # convert to g(C) m^-2
LMA_C[LMA_C>100.]=np.nan

# Create figures
# Figure 1 - vertical distribution of leaf traits:
# Vcmax, Rd, N,
# LMA, Leaf Thickness, C:N ratio
plt.figure(1, facecolor='White',figsize=[9,6])

ax1 = plt.subplot2grid((2,3),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 

mask = ~np.isnan(Vcmax[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Vcmax[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Vcmax[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Vcmax[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])

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

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str
ax1.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax1.set_ylabel('Height / m')
ax1.set_xlabel('$V_{cmax}$')

ax2 = plt.subplot2grid((2,3),(0,1),sharey=ax1)
ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(Rd[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax2.plot(Rd[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Rd[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Rd[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Rd[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Rd[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
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
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_ylabel('Height / m')
ax2.set_xlabel('$R_d$')

ax3 = plt.subplot2grid((2,3),(0,2),sharey=ax1)
ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
#ax3.plot(N[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
#ax3.plot(N[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')
ax3.plot(Narea[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax3.plot(Narea[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

#mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
#m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
#mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
#m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
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
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_ylabel('Height / m')
#ax3.set_xlabel('%N')
ax3.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')


ax4 = plt.subplot2grid((2,3),(1,0),sharey=ax1)
ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax4.plot(LMA_C[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax4.plot(LMA_C[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')
print "LMA(C)"
mask = ~np.isnan(LMA_C[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
print "OG: ", np.mean(LMA_C[ftype=='OG'][mask]), " +/- ", np.std(LMA_C[ftype=='OG'][mask])
m1, c1, r1, p1, err1 = stats.linregress(LMA[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(LMA_C[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
print "SL: ", np.mean(LMA_C[ftype=='SL'][mask]), " +/- ", np.std(LMA_C[ftype=='SL'][mask])
m2, c2, r2, p2, err2 = stats.linregress(LMA[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])

mask = ~np.isnan(LMA_C) & ~np.isnan(LeafHeight)
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

stats_str = 'R$^2$=' + '%.3f' % r1**2 + p1str + '\nR$^2$=' + '%.3f' % r2**2 + p2str
ax4.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax4.set_ylabel('Height / m')
ax4.set_xlabel('LMA(C) / g m$^{-2}$')

ax5 = plt.subplot2grid((2,3),(1,1),sharey=ax1)
ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax5.plot(LeafThickness[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax5.plot(LeafThickness[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(LeafThickness[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(LeafThickness[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(LeafThickness[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(LeafThickness[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
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
ax5.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax5.set_ylabel('Height / m')
ax5.set_xlabel('Leaf Thickness / mm')

ax6 = plt.subplot2grid((2,3),(1,2),sharey=ax1)
ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax6.plot(CNratio[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax6.plot(CNratio[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(CNratio[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(CNratio[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(CNratio[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(CNratio[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
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
ax6.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax6.set_ylabel('Height / m')
ax6.set_xlabel('C:N ratio')

plt.tight_layout()
plt.savefig('vertical_trait_distributions.png')

# Figure 2 - N control on photosynthetic rates and respiration rates
plt.figure(2, facecolor='White',figsize=[9,3])

ax1 = plt.subplot2grid((1,3),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax1.plot(N[ftype=='OG'],Vcmax[ftype=='OG'],'.',color='blue')
ax1.plot(N[ftype=='SL'],Vcmax[ftype=='SL'],'.',color='red')

mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(Vcmax[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],Vcmax[ftype=='OG'][mask])
mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(Vcmax[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],Vcmax[ftype=='SL'][mask])
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

ax1.set_xlabel('%N')
ax1.set_ylabel('$V_{cmax}$')

ax2 = plt.subplot2grid((1,3),(0,1),sharex=ax1)
ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(N[ftype=='OG'],Rd[ftype=='OG'],'.',color='blue')
ax2.plot(N[ftype=='SL'],Rd[ftype=='SL'],'.',color='red')
mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(Rd[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],Rd[ftype=='OG'][mask])
mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(Rd[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],Rd[ftype=='SL'][mask])
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
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_xlabel('%N')
ax2.set_ylabel('$R_d$')

ax3 = plt.subplot2grid((1,3),(0,2),sharex=ax1)
ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax3.plot(N[ftype=='OG'],Vcmax[ftype=='OG']/Rd[ftype=='OG'],'.',color='blue')
ax3.plot(N[ftype=='SL'],Vcmax[ftype=='SL']/Rd[ftype=='SL'],'.',color='red')

mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(Vcmax[ftype=='OG']/Rd[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],(Vcmax[ftype=='OG']/Rd[ftype=='OG'])[mask])
mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(Vcmax[ftype=='SL']/Rd[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],(Vcmax[ftype=='SL']/Rd[ftype=='SL'])[mask])
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
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_xlabel('%N')
ax3.set_ylabel('$V_{cmax}/R_d$')

plt.tight_layout()
plt.savefig('leaf_N_and_photosynthetic_rates.png')

# Figure 3 - N/area control on photosynthetic rates and respiration rates
plt.figure(3, facecolor='White',figsize=[9,3])

ax1 = plt.subplot2grid((1,3),(0,0))
ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax1.plot(Narea[ftype=='OG'],Vcmax[ftype=='OG'],'.',color='blue')
ax1.plot(Narea[ftype=='SL'],Vcmax[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(Vcmax[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],Vcmax[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(Vcmax[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],Vcmax[ftype=='SL'][mask])
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

ax1.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')
ax1.set_ylabel('$V_{cmax}$')

ax2 = plt.subplot2grid((1,3),(0,1),sharex=ax1)
ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax2.plot(Narea[ftype=='OG'],Rd[ftype=='OG'],'.',color='blue')
ax2.plot(Narea[ftype=='SL'],Rd[ftype=='SL'],'.',color='red')
mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(Rd[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],Rd[ftype=='OG'][mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(Rd[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],Rd[ftype=='SL'][mask])
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
ax2.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax2.set_xlabel('[N]$_{area}$ / g(N)m$^{-2}$')
ax2.set_ylabel('$R_d$')

ax3 = plt.subplot2grid((1,3),(0,2),sharex=ax1)
ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax3.plot(Narea[ftype=='OG'],Vcmax[ftype=='OG']/Rd[ftype=='OG'],'.',color='blue')
ax3.plot(Narea[ftype=='SL'],Vcmax[ftype=='SL']/Rd[ftype=='SL'],'.',color='red')

mask = ~np.isnan(Narea[ftype=='OG']) & ~np.isnan(Vcmax[ftype=='OG']/Rd[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(Narea[ftype=='OG'][mask],(Vcmax[ftype=='OG']/Rd[ftype=='OG'])[mask])
mask = ~np.isnan(Narea[ftype=='SL']) & ~np.isnan(Vcmax[ftype=='SL']/Rd[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(Narea[ftype=='SL'][mask],(Vcmax[ftype=='SL']/Rd[ftype=='SL'])[mask])
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
ax3.annotate(stats_str, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']) 

ax3.set_xlabel('[N]$_{area}$ / g(N)cm$^{-2}$')
ax3.set_ylabel('$V_{cmax}/R_d$')

plt.tight_layout()
plt.savefig('leaf_Narea_and_photosynthetic_rates.png')
plt.show()



#=============================================
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

#### Currently not indexing across into both plot and tree tag - need to do both as it seems that in census, tree tags non unique!

subplot = np.zeros(N_branches)*np.nan

plots = np.unique(plot)
for i in range(0,N_branches):
    tree_index = np.all((census_plot==plot[i], np.any((tree_tag == tree[i], alt_tag == tree[i]),axis=0)),axis=0)
    if tree_index.sum() == 1:
        #print "found tree :-)"
        subplot[i] = census_subplot[tree_index]
    elif tree_index.sum()==0:
        print "cannot find tree from traits record in the census data"
        print "\t tree: ", tree[i], "; plot: ", plot[i]
        subplot[i] = np.nan
    else:
        print "we have a problem - more than one tree in this plot has the same tag "
        subplot[i] = np.nan

# load in light transmittance profiles
light_file = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/src/output/BALI_subplot_lighttransmittance.npz'
I = np.load(light_file)

light_availability = np.zeros(N_branches)*np.nan
#currently I has vertical resolution of 1 m.  Bin 0 -> 0-1 m.  Therefore take floor of leaf height to find index required
for i in range(0,N_branches):
    if np.all((np.isfinite(LeafHeight[i]), np.isfinite(subplot[i]))):
        ht_index = np.floor(LeafHeight[i])
        light_availability[i]=I[plot[i]][subplot[i-1],ht_index]
    else:
        light_availability[i]=np.nan                               
        
branch, spp, genus, N, Narea, C, CNratio, SLA, LMA, LeafArea, LeafThickness, LeafHeight, VPD, Rd, Vcmax, Jmax, ShadeTag, ftype

# write traits data to a csv file for ingestion into R
#plot,subplot,forest_type,branch,shade_tag,spp,genus,leaf_height,light_availability,leaf_thickness,leaf_area,LMA,C%,Carea,N%,Narea,Vcmax,Rd
out = open('BALI_leaf_traits.csv','w')
out.write('plot,subplot,forest_type,branch,shade_tag,spp,genus,leaf_height,light_availability,leaf_thickness,leaf_area,LMA,C%,Carea,N%,Narea,Vcmax,Rd\n')
for i in range(0,N_branches):
    out.write(plot[i] + ',' + str(subplot[i])  + ', ' + ftype[i] + ',' + branch[i]  + ',' + str(shade_tag[i])  + ',' + spp[i] + ',' + genus[i] + ',' + str(LeafHeight[i]) + ',' + str(light_availability[i]) + ',' +str( LeafThickness[i]) + ',' +str( LeafArea[i]) + ',' + str(LMA[i]) + ',' + str(C[i]) + str(LMA_C[i]) + ',' + str(N[i]) + ',' + str(Narea[i]) + ',' + str(Vcmax[i]) + ',' + str(Rd[i]) + '\n')

    out.close()


