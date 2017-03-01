# This library contains code to download MODIS data and prepare it for assimilation into CARDAMOM
import numpy as np
import requests

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

    return LAI, LAI_std, LAI_qc
