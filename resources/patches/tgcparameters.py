#!/usr/bin/env python

import sys
import os
import json
import math
from glob import glob
import xarray as xr
import numpy as np

input = sys.argv[1]
key = os.path.basename(input)[:-len(".SAFE")]
ecmwf = glob(os.path.join(input, "GRANULE", "*", "AUX_DATA", "AUX_ECMWFT"))[0]
cams  = glob(os.path.join(input, "GRANULE", "*", "AUX_DATA", "AUX_CAMS??"))[0]
ds1 = xr.open_dataset(ecmwf, engine='cfgrib')
ds2 = xr.open_dataset(cams, engine='cfgrib')
aot = ds2['aod550'].values
u = ds1['u10'].values
v = ds1['v10'].values
lat = ds1['latitude'].values
lon = ds1['longitude'].values
result = {
  'wind': np.sqrt(u*u+v*v).reshape(9*9).tolist(),
  'lon': np.tile(lon, (9,1)).reshape(9*9).tolist(),
  'lat': np.tile(lat, (9,1)).T.reshape(9*9).tolist(),
  'aot': aot.reshape(9*9).tolist(),
}
with open(f"{key}-TGC-parameters.json", "w") as output:
    json.dump(result, output)
print(f"{key}-TGC-parameters.json")
