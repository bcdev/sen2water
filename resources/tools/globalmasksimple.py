#!/usr/bin/env python

print("loading software ... ")

import rioxarray as rio
import xarray as xr
import numpy as np
import scipy
from pyproj import CRS, Transformer
from pyproj.enums import TransformDirection
import math
import sys

if len(sys.argv) != 7:
    print('Usage:')
    print('  globalmasksimple.py <worldcoverpath> <globalislandspath> <continentalshorelinepath> <countslist> <zerogranulepath> <outputpath>')
    print('e.g.')
    print("  globalmasksimple.py /data/worldcover/worldcoverlandmask.zarr /data/globalislands/globalislandsinlandmask.zarr /data/continentalshoreline/continentalshorelinemask.zarr globalmaskcounts.list zero-masks/32UME-zeros.tif global-masks/s2w-globalmask-32UME.tif")
    sys.exit(1)

"""
worldcover_path = "/data/worldcover/worldcoverlandmask.zarr"
globalislands_path = "/data/globalislands/globalislandsinlandmask.zarr"
continentalshoreline_path = "/data/continentalshoreline/continentalshorelinemask.zarr"
zeromask_path = "zero-masks/01CDJ-zeros.tif"
s2w_mask_path = "global-masks/s2w-globalmask-01CDJ.tif"
countslist_path = "/data/sen2water/globalmaskcounts.list"
zeromask_path = "zero-masks/32UME-zeros.tif"
s2w_mask_path = "global-masks/s2w-globalmask-32UME.tif"
"""

worldcover_path = sys.argv[1]
globalislands_path = sys.argv[2]
continentalshoreline_path = sys.argv[3]
countslist_path = sys.argv[4]
zeromask_path = sys.argv[5]
s2w_mask_path = sys.argv[6]

# define some constants and functions

INVALID = 0
OCEAN = 64
INLANDWATER = 96
BOTTOMREFLECTION = 128
COASTAL = 160
LAND = 192

RIVERWATER=3
SHORE_TRANSITION = 2
SHORE_INLANDWATER = 1

FOUR_NEIGHBOURS=np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=np.int8)
EIGHT_NEIGHBOURS=np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=np.int8)

# replace neighbours of value source that are value dest by value insert
def dilate(mask, source_value, neighbours, steps, dest_value, insert_value):
    source_mask = (mask == source_value) | (mask == insert_value)
    dest_mask = mask == dest_value
    dilated = scipy.ndimage.binary_dilation(source_mask, neighbours, steps, dest_mask)
    mask[dest_mask & dilated] = insert_value
    del source_mask, dest_mask, dilated

# replace neighbours of value source that are value dest by value insert
def dilate_range(mask, source_value1, source_value9, neighbours, steps, dest_value, insert_value):
    source_mask = (mask>=source_value1) & (mask<=source_value9)
    dest_mask = mask == dest_value
    dilated = scipy.ndimage.binary_dilation(source_mask, neighbours, steps, dest_mask)
    mask[dest_mask & dilated] = insert_value
    del source_mask, dest_mask, dilated

print("determine bounding box of granule in geographic worldcover grid ...")

hr = rio.open_rasterio(zeromask_path)
hr_width = len(hr.x)
hr_height = len(hr.y)
hr_size = hr_height * hr_width
hr_step = 60.0
hr_left = hr.x.data[0]
hr_top = hr.y.data[0]

wc = xr.open_dataset(worldcover_path, engine='zarr', mask_and_scale=False)

# create geo-trafo WGS84 -> UTM
hr_crs = CRS.from_cf(hr.spatial_ref.attrs)
wc_crs = CRS.from_cf(wc.spatial_ref.attrs)
wgs_to_utm = Transformer.from_crs(wc_crs, hr_crs)

# map border of granule to geographic grid
half_a_pixel = hr_step / 2
# extract image border (top, bottom, left, right) for transformation
hr_y = hr.y.values
hr_x = hr.x.values
hr_y_border = np.stack((np.tile(hr_y[0] + half_a_pixel, (hr_width)),
                        np.tile(hr_y[-1] - half_a_pixel, (hr_width)),
                        hr_y,
                        hr_y))
hr_x_border = np.stack((hr_x,
                        hr_x,
                        np.tile(hr_x[0] - half_a_pixel, (hr_width)),
                        np.tile(hr_x[-1] + half_a_pixel, (hr_width))))
# move four pixels of border into image corners
hr_x_border[0,0] -= half_a_pixel
hr_x_border[0,-1] += half_a_pixel
hr_x_border[1,0] -= half_a_pixel
hr_x_border[1,-1] += half_a_pixel
# transform border from UTM to WGS84 geographic coordinates
hr_lat,hr_lon = wgs_to_utm.transform(hr_x_border.flatten(), hr_y_border.flatten(), direction=TransformDirection.INVERSE)
del hr_y_border, hr_x_border, hr_y, hr_x

# determine bounding box in geographic grid
hr_lat_min = np.min(hr_lat)
hr_lat_max = np.max(hr_lat)
# turn lon by first lon as reference to determine minmax lon of border
hr_lon_ref = hr_lon[0]
hr_lon_min = (np.min((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
hr_lon_max = (np.max((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
del hr_lat, hr_lon
# determine minmax pixel coordinates of worldcover subset
wc_j_min = int(math.floor((90.0 - hr_lat_max) / 180.0 * wc.y.shape[0]))
wc_j_max = int(math.ceil((90.0 - hr_lat_min) / 180.0 * wc.y.shape[0]))
wc_i_min = int(math.floor((hr_lon_min + 180.0) / 360.0 * wc.x.shape[0]))
wc_i_max = int(math.ceil((hr_lon_max + 180.0) / 360.0 * wc.x.shape[0]))
i_step = int(1.0 / math.cos((hr_lat_min + hr_lat_max) / 2 * math.pi / 180.0))
print(f"extent {wc_j_min=} {wc_i_min=} {wc_j_max-wc_j_min=} {wc_i_max-wc_i_min=}")
print(f"latitude subsampling factor for geographic grids: {i_step}")
if hr_lon_min > hr_lon_max:
    print('crossing antimeridian')

print("determine mapping from geographic to UTM grid ...")

# create subset of worldcover land mask
if wc_i_min <= wc_i_max:
    wc_subset = wc.landmask[wc_j_min:wc_j_max, wc_i_min:wc_i_max:i_step]
else:
    wc_subset = xr.concat((wc.landmask[wc_j_min:wc_j_max, wc_i_min:], wc.landmask[wc_j_min:wc_j_max, :wc_i_max]), dim='x')
    wc_subset = wc_subset[:,::i_step]
del wc

# transform land mask coordinates from WGS84 to UTM
wc_lat = np.tile(wc_subset.y.values, (len(wc_subset.x), 1)).T
wc_lon = np.tile(wc_subset.x.values, (len(wc_subset.y), 1))
wc_x, wc_y = wgs_to_utm.transform(wc_lat, wc_lon, direction=TransformDirection.FORWARD)
del wc_lat, wc_lon
# determine pixel coordinates in HR-OC UTM grid for landmask subset
# label pixels outside of HR-OC mask granule
wc_j = ((hr_top - wc_y + hr_step / 2) / hr_step).astype(np.int32)
wc_i = ((wc_x - hr_left + hr_step / 2) / hr_step).astype(np.int32)
del wc_x, wc_y

# map the landmask to the UTM grid
# bin_indices is the flattened wc coordinates mapped to hr grid bin indices
# landmask_values is the corresponding flattened wc values, 0 or 1
is_inside_granule = ((wc_j >= 0) & (wc_j < hr_height) & (wc_i >= 0) & (wc_i < hr_width)).astype(bool)
inside_outside_bin_indices = wc_j * hr_width + wc_i
bin_indices = inside_outside_bin_indices[is_inside_granule]
del inside_outside_bin_indices
all_bin_counts = np.bincount(bin_indices, minlength=hr_size)
del wc_j, wc_i

print("subset and reproject worldcover land mask ...")

wc_values = wc_subset.data[is_inside_granule]
del wc_subset
wc_bin_counts = np.bincount(bin_indices[wc_values!=0], minlength=hr_size)
wc_bin_mask = 2 * wc_bin_counts >= all_bin_counts
del wc_bin_counts
wc_landmask = wc_bin_mask.reshape(hr_height, hr_width).astype(np.uint8)
del wc_bin_mask
"""
wc_da = xr.DataArray(wc_landmask.reshape((1, hr_height, hr_width)), coords=hr.coords, dims=hr.dims, attrs=hr.attrs)
wc_da.rio.to_raster('worldcover-32UME.tif')
"""

print("subset and reproject globalislands land mask ...")

gi = xr.open_dataset(globalislands_path, engine='zarr', mask_and_scale=False)
if wc_i_min <= wc_i_max:
    gi_subset = gi.inlandmask[wc_j_min:wc_j_max, wc_i_min:wc_i_max:i_step]
else:
    gi_subset = xr.concat((gi.inlandmask[wc_j_min:wc_j_max, wc_i_min:], gi.inlandmask[wc_j_min:wc_j_max, :wc_i_max]), dim='x')
    gi_subset = gi_subset[:,::i_step]
del gi
gi_values = gi_subset.data[is_inside_granule]
del gi_subset
gi_bin_counts = np.bincount(bin_indices[gi_values!=0], minlength=hr_size)
gi_bin_mask = 2 * gi_bin_counts >= all_bin_counts
del gi_bin_counts
gi_landmask = gi_bin_mask.reshape(hr_height, hr_width).astype(np.uint8)
del gi_bin_mask

print("subset and reproject continentalshoreline land mask ...")

cs = xr.open_dataset(continentalshoreline_path, engine='zarr', mask_and_scale=False)
if wc_i_min <= wc_i_max:
    cs_subset = cs.continentalmask[wc_j_min:wc_j_max, wc_i_min:wc_i_max:i_step]
else:
    cs_subset = xr.concat((cs.continentalmask[wc_j_min:wc_j_max, wc_i_min:], cs.continentalmask[wc_j_min:wc_j_max, :wc_i_max]), dim='x')
    cs_subset = cs_subset[:,::i_step]
del cs
cs_values = cs_subset.data[is_inside_granule]
del cs_subset
del is_inside_granule
cs_bin_counts = np.bincount(bin_indices[cs_values!=0], minlength=hr_size)
cs_bin_mask = 2 * cs_bin_counts >= all_bin_counts
del cs_bin_counts, all_bin_counts
cs_landmask = cs_bin_mask.reshape(hr_height, hr_width).astype(np.uint8)
del cs_bin_mask

print("combine masks ...")

# odd=wc land, 2367=gi land, 4 upwards=cs land, i.e. 0=OCEAN, odd1357=LAND, 26=INLANDWATER, 4=RIVERWATER
"""
mask = wc_landmask + 2 * gi_landmask + 4 * cs_landmask
combinedmask_da = xr.DataArray(mask.reshape((1, hr_height, hr_width)), coords=hr.coords, dims=hr.dims, attrs=hr.attrs)
combinedmask_da.rio.to_raster('combined-01CDJ.tif')
"""

mask = wc_landmask
mask[(mask != 0)] = LAND
mask[((mask == 0) & (gi_landmask!=0))] = INLANDWATER
mask[((mask == 0) & (cs_landmask==0))] = OCEAN
mask[(mask==0)] = RIVERWATER
del wc_landmask, gi_landmask, cs_landmask

print("buffering ocean into rivers and back to preserve ocean included in continental shoreline ...")

depth = 24  # 32
for i in range(depth):
    dilate(mask, OCEAN, FOUR_NEIGHBOURS, 1, RIVERWATER, OCEAN)

# ... and by reverting in river mouths
for i in range(depth):
    dilate(mask, RIVERWATER, FOUR_NEIGHBOURS, 1, OCEAN, RIVERWATER)

mask[mask==RIVERWATER] = INLANDWATER

print("buffering ocean into inlandwater and back to extend ocean up to worldcover coastline ...")

# Are there lakes close to shoreline that are partially marked as ocean?
# buffer ocean into inlandwater by n pixels
# buffer inlandwater back into ocean in estuaries
depth = 48  # 24
for i in range(depth):
    dilate(mask, OCEAN, FOUR_NEIGHBOURS, 1, INLANDWATER, OCEAN)

# ... and by reverting in river mouths
for i in range(depth):
    dilate(mask, INLANDWATER, FOUR_NEIGHBOURS, 1, OCEAN, INLANDWATER)

print("buffering ocean into land to mark coastal land ...")

depth = 28
for i in range(depth):
    dilate(mask, OCEAN, EIGHT_NEIGHBOURS, 1, LAND, COASTAL)

# dilate from transition zone into coastal and from inland water into coastal
# temporarily label shore of both

print("mask transition zone from ocean to inland water ...")

# values 64 < v < 96 are the transition zone from ocean to inland water
# incremental values 7+distance from ocean up to 32 pixels
# 4 neighbours for the first two steps to avoid jumps over locks
depth = 32
for i in range(depth):
    dilate(mask, OCEAN+i, FOUR_NEIGHBOURS if i<2 else EIGHT_NEIGHBOURS, 1, INLANDWATER, OCEAN+i+1)

shore_width = 4
transition_depth = 32
for i in range(shore_width):
    if i == 0:
        dilate_range(mask, OCEAN+1, INLANDWATER-1, EIGHT_NEIGHBOURS, 1, COASTAL, SHORE_TRANSITION)
    else:
        dilate(mask, SHORE_TRANSITION, EIGHT_NEIGHBOURS, 1, COASTAL, SHORE_TRANSITION)
    dilate(mask, INLANDWATER, EIGHT_NEIGHBOURS, 1, COASTAL, SHORE_INLANDWATER)

# count distance from ocean in shore of transition zone
for i in range(transition_depth):
    dilate(mask, OCEAN if i == 0 else COASTAL+i, EIGHT_NEIGHBOURS, 1, SHORE_TRANSITION, COASTAL+i+1)

# revert inland water shore and unreached transition zone shore
mask[(mask==SHORE_INLANDWATER)|(mask==SHORE_TRANSITION)] = COASTAL

print("counting pixel classes ...")

land_count = np.count_nonzero(mask>=128)
inlandwater_count = np.count_nonzero((mask<128) & (mask > 64))
ocean_count = np.count_nonzero(mask<=64)
with open(countslist_path, "a") as f:
    f.write(f"{s2w_mask_path}\t{ocean_count}\t{inlandwater_count}\t{land_count}\n")

print("write s2w mask ...")

# write result to s2w mask
s2w = xr.DataArray(mask.reshape((1, hr_height, hr_width)), coords=hr.coords, dims=hr.dims, attrs=hr.attrs)
s2w.rio.to_raster(s2w_mask_path, compress='LZW', tiled=True)

print(s2w_mask_path)
