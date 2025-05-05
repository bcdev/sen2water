#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Converts GlobalIslands into a zarr"""

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

zarrdir = '/data/globalislands/globalislandsinlandmask.zarr'

d2 = 360.0 / (36000 * 120 * 2)

path = sys.argv[1]
name = path[path.rfind('/')+1:]  # globalislands-h007v07.tif

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
        ds = xr.Dataset({'inlandmask': da, 'spatial_ref': spatial_ref},
                        attrs={'title': 'WorldCover 10m land mask'})
        ds.to_zarr(zarrdir, mode='w', compute=False)
        print(zarrdir, 'initialised')

# globalislands-h007v07.tif
southing = int(name[19:21]) * 3
west_easting = int(name[15:18]) * 3
upper_left_y = southing * 12000
upper_left_x = west_easting * 12000

with xr.open_dataset(path) as ds_lc:
    mask = ds_lc.band_data.values[0]
    da_mask = xr.DataArray(mask, coords={'y': ds_lc.y, 'x': ds_lc.x}, dims=['y', 'x'], attrs={'fill_value': 0})
    ds_mask = xr.Dataset({'inlandmask': da_mask}).drop_vars(['y', 'x'])
    ds_mask_slice = ds_mask.isel({'y':slice(0,36000),'x':slice(0,36000)})
    ds_mask_slice.to_zarr(zarrdir, mode='r+', region={'y': slice(upper_left_y, upper_left_y+36000), 'x': slice(upper_left_x, upper_left_x+36000)})

print(upper_left_y // 3 // 12000, upper_left_x // 3 // 12000, 'filled')
