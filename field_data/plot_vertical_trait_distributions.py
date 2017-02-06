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



chem_file = '../../../BALI_traits_data/CombinedPlots/BALI_traits_CN_02122016.csv'
spp_file = '../../../BALI_traits_data/CombinedPlots/BALI_species_02122016.csv'
photo_file = '../../../BALI_traits_data/CombinedPlots/Photosynthesis_Combined_14122016.csv'
branch_file = '../../../BALI_traits_data/CombinedPlots/ParameterTrees_Combined_14122016.csv'
leaf_file = '../../../BALI_traits_data/CombinedPlots/LeafArea_Combined_14122016.csv'

# Load in the traits
branch, spp, genus, N, C, CNratio, SLA, LMA, LeafArea, LeafThickness, LeafHeight, VPD, Rd, Vcmax, Jmax, ShadeTag, ftype = field.collate_branch_level_traits(chem_file,photo_file,branch_file,leaf_file,spp_file)

LeafThickness[LeafThickness>0.60]=np.nan
Vcmax[Vcmax>150]=np.nan
LMA[LMA>250]=np.nan
N[N>4]=np.nan



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
ax3.plot(N[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax3.plot(N[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(N[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(N[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(N[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(N[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])
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
ax3.set_xlabel('%N')

ax4 = plt.subplot2grid((2,3),(1,0),sharey=ax1)
ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2) 
ax4.plot(LMA[ftype=='OG'],LeafHeight[ftype=='OG'],'.',color='blue')
ax4.plot(LMA[ftype=='SL'],LeafHeight[ftype=='SL'],'.',color='red')

mask = ~np.isnan(LMA[ftype=='OG']) & ~np.isnan(LeafHeight[ftype=='OG'])
m1, c1, r1, p1, err1 = stats.linregress(LMA[ftype=='OG'][mask],LeafHeight[ftype=='OG'][mask])
mask = ~np.isnan(LMA[ftype=='SL']) & ~np.isnan(LeafHeight[ftype=='SL'])
m2, c2, r2, p2, err2 = stats.linregress(LMA[ftype=='SL'][mask],LeafHeight[ftype=='SL'][mask])

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
ax4.set_xlabel('LMA')

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
plt.show()
