import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

def downsample_with_randomforest(rs_variables,target_variable):
    print "\t downsampling with random forest regression"
    # first of all get rid of nodata values in target data to build random forest
    rs = rs_variables[:,np.isfinite(target_variable)]
    target = target_variable[np.isfinite(target_variable)]

    # split the data into calibration and validation datasets
    rs_cal, rs_val, target_cal, target_val = train_test_split(rs,target,train_size=0.5)

    # build random forest regressor
    n_trees_in_forest = 50
    min_samples_leaf = 50.  # increase this if validation quality much worse than calibration
    n_cores = 5
    randomforest = RandomForestRegressor(n_estimators=n_trees_in_forest, min_samples_leaf=min_samples_leaf,  bootstrap=True, oob_score=True, n_jobs=n_cores)
    randomforest.fit(rs_cal,target_cal)

    print "\t\tcalibration score = %.3f" randomforest.score(rs_cal,target_cal)
    print "\t\tvalidation score = %.3f" randomforest.score(rs_val,target_val)
    rs_rf = randomforest.predict(rs_variables)
    return rs_rf
