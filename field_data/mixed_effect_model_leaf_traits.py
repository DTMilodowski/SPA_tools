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
def normalise_column(df_column):
    pd_column = (df_column-df_column.mean())/df_column.std()
    return df_column

# mask nodata variables
def mask_nodata(df,variables):
    data_locs = np.zeros((df.shape[0],len(variables)),dtype = bool)
    for vv in range(0,len(variables)):
        print variables[vv]
        data_locs[:,vv] = np.isfinite(df[variables[vv]])
    nodata_mask = np.all(data_locs,axis=1)   
    return nodata_mask

# run fixed effects model
def run_mixed_effects_model(data,variables_to_use,model,group_vars):
    indices = mask_nodata(data,variables_to_use)
    data_sub = data[indices].copy()

    # scale the data columns that you are going to use so that mean = 0 and standard deviation =1
    for vv in range(0,len(variables_to_use)):
        print variables_to_use[vv]
        data_sub[variables_to_use[vv]]=normalise_column(data_sub[variables_to_use[vv]])


    md = smf.mixedlm(model, data_sub, groups=data_sub[group_vars])
    mdf = md.fit()
    print mdf.summary()
    return md,mdf

# run simple linear model
def run_ols(data,variables_to_use,model):
    indices = mask_nodata(data,variables_to_use)
    data_sub = data[indices].copy()

    # scale the data columns that you are going to use so that mean = 0 and standard deviation =1
    for vv in range(0,len(variables_to_use)):
        print variables_to_use[vv]
        data_sub[variables_to_use[vv]]=normalise_column(data_sub[variables_to_use[vv]])


    md = smf.ols(model, data_sub)
    mdf = md.fit()
    print mdf.summary()
    return md,mdf

#==============================================================
    
# traits file
leaf_traits_file = 'BALI_leaf_traits.csv'
data = pd.read_csv(leaf_traits_file)

# extract subset of data containing data values for variables of interest 
variables_to_use  = ['Vcmax', 'Narea','Carea','leaf_height','light_availability','Rd']
variables_to_use  = ['Narea','leaf_height','light_availability']

# choose model
model = "Narea ~ leaf_height+light_availability"
groups = "spp"

# run model
md, mdf = run_mixed_effects_model(data,variables_to_use,model,groups)








indices = mask_nodata(data,variables_to_use)
data_sub = data[indices].copy()








# scale the data columns that you are going to use so that mean = 0 and standard deviation =1
for vv in range(0,len(variables_to_use)):
    print variables_to_use[vv]
    data_sub[variables_to_use[vv]]=normalise_column(data_sub[variables_to_use[vv]])


md = smf.mixedlm("Rd ~ Narea+leaf_height+light_availability", data[indices], groups=data[indices]["spp"])
mdf = md.fit()
mdf.summary()

run_mixed_effects_model(data,variables_to_use,model,groups=data[indices]["spp"])

# simple linear model, no categoricals
md = smf.ols("Rd ~ Narea+leaf_height+light_availability", data[indices])
mdf = md.fit()

