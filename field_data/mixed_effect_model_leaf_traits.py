import numpy as np
from matplotlib import pyplot as plt
import sys

import statsmodels.api as sm
import statsmodels.formula.api as smf

import pandas as pd


from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2


# define a function to normalise a dataset so that mean in zero and standard deviation is one. 
# This allows easier comparison of effect size
def normalise_column(pd_column):
    pd_column = (pd_column-pd_column.mean())/pd_column.std()
    return pd_column

# traits file
leaf_traits_file = 'BALI_leaf_traits.csv'

#datatype={'names':('plot','subplot','forest_type','branch','shade_tag','spp','genus','leaf_height','light_availability','leaf_thickness','leaf_area','LMA','C%','Carea','N%','Narea','Vcmax','Rd'),'formats':('S16','i8','S16','S16','i8','S32','S32','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16')}
#data = np.genfromtxt(chem_file,skiprows=1,delimiter=',',dtype=datatype)

data = pd.read_csv(leaf_traits_file)

indices = np.all((np.isfinite(data['Vcmax']),np.isfinite(data['Narea'])),axis =0)

# extract subset of data containing data values for variables of interest 
data_sub = data[indices].copy()

# scale the data columns that you are going to use so that mean = 0 and standard deviation =1
variables_to_use  = ['Vcmax', 'Narea','Carea','leaf_height','light_availability','Rd']
for vv in range(0,len(variables_to_use)):
    print variables_to_use[vv]
    data_sub[variables_to_use[vv]]=normalise_column(data_sub[variables_to_use[vv]])


md = smf.mixedlm("Vcmax ~ Narea", data[indices], groups=data[indices]["spp"])
mdf = md.fit()


