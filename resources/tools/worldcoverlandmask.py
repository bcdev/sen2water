# -*- coding: utf-8 -*-

"""
Tool to convert WorldCover into a zarr

Usage:
ssh cate
cd /home/boe
. /opt/miniconda-sen2water-0.2/bin/activate
export PATH=$(pwd):/opt/calvalus-2.3/bin:$PATH
export PYTHONPATH=/opt/calvalus-2.3/lib:/data/sen2water/processor-test/hrocresampling
for f in /data/worldcover/ESA_WorldCover*tif; do
    worldcoverlandmask.py $f
ls -l /data/worldcover/worldcoverlandmask.zarr
"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "proprietary, TBD"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import numpy as np
import dask.array as da
import xarray as xr
import sys
import os

datadir = '/data/worldcover'
zarrdir = datadir + '/worldcoverlandmask.zarr'

d2 = 360.0 / (36000 * 120 * 2)

path = sys.argv[1]
name = path[path.rfind('/')+1:]

if not os.path.exists(zarrdir):
    with xr.open_dataset(path) as ds_lc:
        spatial_ref = ds_lc.spatial_ref
        spatial_ref.attrs['GeoTransform'] = [ -180.0, 2*d2, 0.0, 90.0, 0.0, -2*d2 ]
        x = da.arange(-180.0+d2, 180.0, 2*d2, chunks=36000, dtype=np.float64)
        y = da.arange(90-d2, -90, -2*d2, chunks=36000, dtype=np.float64)
        data = da.zeros((36000*60, 36000*120), chunks=(36000, 36000), dtype=np.byte)
        da = xr.DataArray(data,
                          coords={'y': y, 'x': x}, 
                          dims=['y', 'x'], name='landmask', 
                          attrs={'fill_value': 0, '_FillValue': 0, 'flag_values': [0, 1], 'flag_meanings': ['water', 'land']})
        ds = xr.Dataset({'landmask': da, 'spatial_ref': spatial_ref},
                        attrs={'title': 'WorldCover 10m land mask'})
        ds.to_zarr('worldcoverlandmask.zarr', mode='w', compute=False)
        print(zarrdir, 'initialised')

# ESA_WorldCover_10m_2021_v200_N36E105_Map.tif
northing = int(name[30:32])
easting = int(name[33:36])
upper_left_y = (87 - northing if name[29] == 'N' else 90 + northing) * 12000
upper_left_x = (180 - easting if name[32] == 'W' else 180 + easting) * 12000

with xr.open_dataset(path) as ds_lc:
    lc_class = ds_lc.band_data.values[0]
    mask = ((lc_class != 80) & (np.isnan(lc_class) == 0)).astype(np.int8)
    da_mask = xr.DataArray(mask, coords={'y': ds_lc.y, 'x': ds_lc.x}, dims=['y', 'x'], attrs={'fill_value': 0})
    ds_mask = xr.Dataset({'landmask': da_mask}).drop_vars(['y', 'x'])
    ds_mask_slice = ds_mask.isel({'y':slice(0,36000),'x':slice(0,36000)})
    ds_mask_slice.to_zarr(zarrdir, mode='r+', region={'y': slice(upper_left_y, upper_left_y+36000), 'x': slice(upper_left_x, upper_left_x+36000)})

print(upper_left_y // 3 // 12000, upper_left_x // 3 // 12000, 'filled')
