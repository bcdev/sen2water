#!/usr/bin/env python

print("loading software ... ")

import rioxarray as rio
import xarray as xr
import numpy as np
from pyproj import CRS, Transformer
from pyproj.enums import TransformDirection
import sys

if len(sys.argv) != 5:
    print('Usage:')
    print('  landcovergranule.py <landcoverpath> <countslist> <zerogranulepath> <outputpath>')
    print('e.g.')
    print("  landcovergranule.py /data/fire/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.tif lc-counts.list zero-masks/32UME-zeros.tif lc-masks/lcmask-32UME.tif")
    sys.exit(1)

"""
landcover_path = "/data/fire/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.tif"
countslist_path = "lccounts.list"
#zeromask_path = "zero-masks/01CDJ-zeros.tif"
zeromask_path = "zero-masks/32UME-zeros.tif"
#lc_mask_path = "lc-masks/lc-globalmask-01CDJ.tif"
lc_mask_path = "lc-masks/lc-globalmask-32UME.tif"
"""

landcover_path = sys.argv[1]
countslist_path = sys.argv[2]
zeromask_path = sys.argv[3]
lc_mask_path = sys.argv[4]

print("determine mapping from LC to granule ...")

hr = rio.open_rasterio(zeromask_path)
hr_width = len(hr.x)
hr_height = len(hr.y)
hr_size = hr_height * hr_width
hr_step = 60.0
hr_left = hr.x.data[0]
hr_top = hr.y.data[0]

lc = rio.open_rasterio(landcover_path)

# create geo-trafo WGS84 -> UTM
hr_crs = CRS.from_cf(hr.spatial_ref.attrs)
lc_crs = CRS.from_cf(lc.spatial_ref.attrs)
wgs_to_utm = Transformer.from_crs(lc_crs, hr_crs)

# map granule to geographic grid
hr_y = hr.y.values
hr_x = hr.x.values
hr_y_m = np.tile(hr_y, (hr_x.shape[0],1)).T
hr_x_m = np.tile(hr_x, (hr_y.shape[0],1))
# transform border from UTM to WGS84 geographic coordinates
hr_lat,hr_lon = wgs_to_utm.transform(hr_x_m.flatten(), hr_y_m.flatten(), direction=TransformDirection.INVERSE)
del hr_y_m, hr_x_m, hr_y, hr_x

print("reprojecting LC to granule with nearest neighbour resampling ...")

# convert to nearest neighbour pixel coordinates in lc grid and resample to lc hr grid
lc_height = lc.shape[1]
lc_width = lc.shape[2]
hr_j = ((90.0 - hr_lat) / 180.0 * lc_height).astype(np.int64)
hr_i = ((hr_lon + 180.0) / 360.0 * lc_width).astype(np.int32)
hr_ij = hr_j * lc_width + hr_i
hr_lc = (lc.values.reshape((lc_height * lc_width))[hr_ij]).reshape((hr_height, hr_width))

burnable = ((hr_lc > 0) & (hr_lc <= 180)).astype(np.uint8)

"""
lc_da = xr.DataArray(lc_landmask.reshape((1, hr_height, hr_width)), coords=hr.coords, dims=hr.dims, attrs=hr.attrs)
lc_da.rio.to_raster('worldcover-32UME.tif')
"""

print("counting pixel classes ...")

burnable_count = np.count_nonzero(burnable > 0)
unburnable_count = np.count_nonzero(burnable == 0)
with open(countslist_path, "a") as f:
    f.write(f"{lc_mask_path}\t{unburnable_count}\t{burnable_count}\n")

print("writing burnable mask ...")

# write result to s2w mask
burnable_da = xr.DataArray(burnable.reshape((1, hr_height, hr_width)), coords=hr.coords, dims=hr.dims, attrs=hr.attrs)
burnable_da.rio.to_raster(lc_mask_path, compress='LZW', tiled=True)

print(lc_mask_path)
