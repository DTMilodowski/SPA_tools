# This library contains code to download MODIS data and prepare it for assimilation into CARDAMOM
import numpy as np
import requests

# This function loads in MODIS time series for point locations from the data file downloaded via AppEARS:  https://lpdaacsvc.cr.usgs.gov/appeears/
def load_point_MODIS_LAI_time_series_from_file(MODIS_file):
    dtype={'names':('category','plot','lat','lon','date','MODIS_tile','MOD15A2H_006_Line_Y_500m','MOD15A2H_006_Sample_X_500m','MOD15A2H_006_FparExtra_QC','MOD15A2H_006_FparStdDev_500m','MOD15A2H_006_Fpar_500m','MOD15A2H_006_LaiStdDev_500m','MOD15A2H_006_Lai_500m','MOD15A2H_006_FparLai_QC', 'MOD15A2H_006_FparLai_QC_bitmask','MOD15A2H_006_FparLai_QC_MODLAND','MOD15A2H_006_FparLai_QC_MODLAND_Description','MOD15A2H_006_FparLai_QC_Sensor','MOD15A2H_006_FparLai_QC_Sensor_Description','MOD15A2H_006_FparLai_QC_DeadDetector','MOD15A2H_006_FparLai_QC_DeadDetector_Description', 'MOD15A2H_006_FparLai_QC_CloudState','MOD15A2H_006_FparLai_QC_CloudState_Description','MOD15A2H_006_FparLai_QC_SCF_QC','MOD15A2H_006_FparLai_QC_SCF_QC_Description'),'formats':('S8','S8','f16','f16','S10','S8','i8','i8','i8','f16','f16','f16','f16','i8','S10','S4','S132','S4','S5','S4','S132','S4','S64','S5','S132')}
    data = np.genfromtxt(MODIS_file,skiprows=1,delimiter=',',dtype=dtype)
    return dates, LAI, LAI_std

"""
### Issues with script-based interfacing to LP DAAC server. 

# This is the URL for downloading MODIS data
SERVICES_URL = 'https://lpdaacsvc.cr.usgs.gov/services/appeears-api'


# function to download point time series. start date and end date should be strings using standard date
# format i.e. dd-mm-yyyy
def download_MODIS_LAI_for_point(latitude,longitude,start_date,end_date,ProductAndVersion=u'MOD15A2H.006'):
    
    product = next(p for p in products if p['ProductAndVersion']==ProductAndVersion)
    LAI_layers_req = requests.get('%s/product/%s?format=json' %(SERVICES_URL, product['ProductAndVersion']))
    LAI_layers = LAI_layers_req.json()
    LAI_layer = LAI_layers[u'Lai_500m']
    LAI_std_layer = LAI_layers[u'LaiStdDev_500m']
    LAI_qc_layer = LAI_layers[u'FparLai_QC']

    # Create a variable containing the URL to the sample service
    sample_url = '{0}/sample?'.format(SERVICES_URL)

    # The sample service requires the following arguments: the MODIS product (e.g., 
    # MOD11A1.005), the product layer (e.g., LST_Day_1km), a start and end date, the point
    # location (i.e., latitude and longitude), and the format of the returned data (we only
    # support json at this time).

    sample_args = {
        'product': product['ProductAndVersion'],
        'layer': LAI_layer['Layer'],
        'startdate': start_date,
        'enddate': end_date,
        'coordinate': '{0},{1}'.format(latitude, longitude),
        'format': 'json'
    }
    
    sample_req = requests.get(sample_url, params = sample_args)
    samples = sample_req.json()

    # read LAI and LAI_std into arrays for dates when quality flag is good. 

    return dates, LAI, LAI_std
"""
