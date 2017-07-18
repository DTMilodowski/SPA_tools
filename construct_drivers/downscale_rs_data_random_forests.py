import numpy as np
from sklearn.ensemble import RandomForestRegressor

def downsample_with_randomforest(rs_variables,target_variable):
    
    # first of all get rid of nodata values in target data to build random forest
    rs_train = rs_variables[:,np.isfinite(target_variable)]
    target_train = target_variable[np.isfinite(target_variable)]

    # build random forest regressor
    n_trees_in_forest = 
    max_depth = 
    n_cores = 
    randomforest = RandomForestRegressor(n_estimators=n_trees_in_forest, max_depth=max_depth,  min_weight_fraction_leaf=0.0, max_features='auto', bootstrap=True, oob_score=True, n_jobs=n_cores)
    randomforest.fit(rs_train,target_train)
    print "score = %.3f" randomforest.score(rs_train,target_train)
    rs_rf = randomforest.predict(rs_variables)
