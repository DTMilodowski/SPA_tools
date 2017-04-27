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

# traits file
leaf_traits_file = 'BALI_leaf_traits.csv'

#datatype={'names':('plot','subplot','forest_type','branch','shade_tag','spp','genus','leaf_height','light_availability','leaf_thickness','leaf_area','LMA','C%','Carea','N%','Narea','Vcmax','Rd'),'formats':('S16','i8','S16','S16','i8','S32','S32','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16')}
#data = np.genfromtxt(chem_file,skiprows=1,delimiter=',',dtype=datatype)

data = pd.read_csv(leaf_traits_file)

data_norm = (data - data.mean())/(data.std())

indices = np.all((np.isfinite(data['Vcmax']),np.isfinite(data['Narea'])),axis =0)
md = smf.mixedlm("Vcmax ~ Narea", data[indices], groups=data[indices]["spp"])
mdf = md.fit()


