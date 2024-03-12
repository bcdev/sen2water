#!/data/opt/miniforge3/envs/xcube/bin/python -u

"""
ssh boe@cate
cd /data/sen2water
eval "$(/data/opt/miniforge3/bin/conda shell.bash hook)"
conda activate xcube
s2wmasksimple.py /data/worldcover/worldcoverlandmask.zarr hroc-mask/4nws/raster_mask_tile_32UNF_4nws_32632.tif hroc-mask/3bal/raster_mask_tile_32UNF_3bal_32632.tif
"""

import sys
import math
import rioxarray as rx
import xarray as xr
import scipy
import numpy as np
from pyproj import CRS, Transformer
from pyproj.enums import TransformDirection

"""
Method to convert a HR-OC static mask to a Sen2Water static mask
for all European coastal granules processed by the HR-OC
service. Main purpose of the mask is to distinguish ocean water
from inland water (and land). In the transition zone between ocean
and inland waters at estuaries the distance from the fictive
coastline is encoded into the mask for smooth switching between
different algorithms.

Algorithms used are
 o proj for coordinate transformation between WGS84 and UTM for the worldcover land mask, 
 o scipy's binned_statistic to do the actual resampling in this transformation, and
 o scipy's binary_dilation to fill masks with values, either incrementally or complete.

The main codes are land 192, coastal land 160, inland water 96, ocean 64.
Details:

val binary   zone (colour)
=== ======== =============
192 11000000 land (green)
191 10111111 coastal transition (only near transition water)
... ...      ...
161 10100001 coastal transition (-"-)
160 10100000 coastal land (yellow)
159 10011111 coastal land candidate for coastal inland water (?)
128 10000000 bottom reflection (red)
096 01100000 inland water (lightblue)
095 01011111 transition water
... ...      ...
065 01000001 transition water
064 01000000 ocean (darkblue)
000 00000000 invalid (black)
"""

if len(sys.argv) < 3 or len(sys.argv) > 4:
    print('Usage:')
    print('  s2wmasksimple.py <worldcoverpath> <hrocgranulepath> [<hrocgranulepath>]')
    print('e.g.')
    print("  s2wmasksimple.py /data/worldcover/worldcoverlandmask.zarr hroc-mask/4nws/raster_mask_tile_32UNF_4nws_32632.tif hroc-mask/3bal/raster_mask_tile_32UNF_3bal_32632.tif")
    sys.exit(1)

"""
worldcover_path = '/data/worldcover/worldcoverlandmask.zarr'
hroc_mask_path = '/data/sen2water/hroc-mask/4nws/raster_mask_tile_32UMF_4nws_32632.tif'
hroc_mask2_path = None
"""
worldcover_path = sys.argv[1]
hroc_mask_path = sys.argv[2]
hroc_mask2_path = sys.argv[3] if len(sys.argv) > 3 else None

# define some constants and functions

INVALID = 0
OCEAN = 64
INLANDWATER = 96
BOTTOMREFLECTION = 128
COASTAL = 160
LAND = 192

four_neighbours=np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=np.int8)
eight_neighbours=np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=np.int8)

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

# majority(np.array([0,0,1,1,1])) = 1
# majority(np.array([0,0,0,1,1])) = 0
# majority(np.array([])) = 0
def majority(a):
    return (np.sum(a) + len(a) // 2) // len(a) if len(a) > 0 else 0

# ----------------------------------------

print("prepare HR-OC mask ...")

# read HR-OC mask into mask that we will update and finally write as s2w mask
filename = hroc_mask_path[hroc_mask_path.rfind('/')+1:]
granule = filename[len('raster_mask_tile_'):len('raster_mask_tile_32UME')]
dataset = rx.open_rasterio(hroc_mask_path)
mask = dataset.values[0].astype(np.uint8)
hr_width = len(dataset.x)
hr_height = len(dataset.y)
hr_size = hr_height * hr_width
hr_step = 60.0
hr_left = dataset.x.data[0]
hr_top = dataset.y.data[0]

if hroc_mask2_path:
    print("combine second HR-OC mask ...")
    dataset2 = rx.open_rasterio(hroc_mask2_path)
    mask2 = dataset2.values[0].astype(np.uint8)
    # fill areas empty in mask with values from mask2
    mask[(mask==0)&(mask2!=0)] = mask2[(mask==0)&(mask2!=0)]
    # overwrite land in mask with coastal in mask2
    mask[(mask==1)&(mask2==3)] = 3
    del mask2, dataset2

# ----------------------------------------

print("fill land and ocean in mask ...")

mask[mask==1] = LAND
mask[mask==2] = OCEAN
mask[mask==3] = COASTAL
mask[mask==4] = BOTTOMREFLECTION

# repair invalid border of up to four pixels
# dilate land (192), ocean (64), coastal (160), bottom reflection (128) into invalid (0)
steps = 4
for i in range(steps):
    for v in [LAND, OCEAN, OCEAN, COASTAL, BOTTOMREFLECTION]:
        dilate(mask, v, four_neighbours, 1, INVALID, v)
# fill all invalid areas (0) with land (192) or ocean (64)
for v in [LAND, OCEAN]:
    dilate(mask, v, eight_neighbours, 0, INVALID, v)

# ----------------------------------------

print("determine granule extent and subset worldcover land mask ...")

# read worldcover land mask
wc = xr.open_dataset(worldcover_path, engine='zarr', mask_and_scale=False)

# create geo-trafo WGS84 -> UTM
hr_crs = CRS.from_cf(dataset.spatial_ref.attrs)
wc_crs = CRS.from_cf(wc.spatial_ref.attrs)
wgs_to_utm = Transformer.from_crs(wc_crs, hr_crs)

# transform 2D coordinate arrays of granule from UTM to WGS84
hr_y = np.tile(dataset.y.values, (hr_width, 1))
hr_x = np.tile(dataset.x.values, (hr_height, 1)).T
hr_lat,hr_lon = wgs_to_utm.transform(hr_x, hr_y, direction=TransformDirection.INVERSE)
del hr_y, hr_x

# determine minmax latlon of border
hr_lat_min = min(np.min(hr_lat[0]),
                 np.min(hr_lat[-1]),
                 np.min(hr_lat[:,0]),
                 np.min(hr_lat[:,-1]))
hr_lat_max = max(np.max(hr_lat[0]),
                 np.max(hr_lat[-1]),
                 np.max(hr_lat[:,0]),
                 np.max(hr_lat[:,-1]))
# turn lon by first lon as reference to determine minmax latlon
hr_lon_ref = hr_lon[0,0]
hr_lon_min = (min(np.min((hr_lon[0] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.min((hr_lon[-1] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.min((hr_lon[:,0] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.min((hr_lon[:,-1] - hr_lon_ref + 180.0) % 360.0 - 180.0)) + hr_lon_ref + 180.0) % 360.0 - 180.0
hr_lon_max = (max(np.max((hr_lon[0] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.max((hr_lon[-1] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.max((hr_lon[:,0] - hr_lon_ref + 180.0) % 360.0 - 180.0),
                  np.max((hr_lon[:,-1] - hr_lon_ref + 180.0) % 360.0 - 180.0)) + hr_lon_ref + 180.0) % 360.0 - 180.0

# determine minmax pixel coordinates of worldcover subset
wc_j_min = int(math.floor((90.0 - hr_lat_max) / 180.0 * wc.y.shape[0]))
wc_j_max = int(math.ceil((90.0 - hr_lat_min) / 180.0 * wc.y.shape[0]))

if hr_lon_min > hr_lon_max:
    print('crossing antimeridian - not yet implemented')
    sys.exit(1)

wc_i_min = int(math.floor((hr_lon_min + 180.0) / 360.0 * wc.x.shape[0]))
wc_i_max = int(math.ceil((hr_lon_max + 180.0) / 360.0 * wc.x.shape[0]))

# create subset of worldcover land mask
landmask_wgs84 = wc.landmask[wc_j_min:wc_j_max,wc_i_min:wc_i_max]
"""
landmask_wgs84.rio.to_raster('worldcover-subset.tif')
"""
del hr_lat, hr_lon, wc

print("collocate worldcover land mask ... ... ...")

# transform land mask coordinates from WGS84 to UTM
wc_lat = np.tile(landmask_wgs84.y.values, (len(landmask_wgs84.x), 1)).T
wc_lon = np.tile(landmask_wgs84.x.values, (len(landmask_wgs84.y), 1))
wc_x, wc_y = wgs_to_utm.transform(wc_lat, wc_lon, direction=TransformDirection.FORWARD)
del wc_lat, wc_lon

# determine pixel coordinates in HR-OC UTM grid for landmask subset
# label pixels outside of HR-OC mask granule
wc_j = ((hr_top - wc_y + hr_step / 2) / hr_step).astype(np.int32)
wc_i = ((wc_x - hr_left + hr_step / 2) / hr_step).astype(np.int32)
del wc_x, wc_y
wc_j[(wc_j < 0)|(wc_j >= hr_height)|(wc_i < 0)|(wc_i >= hr_width)] = -1
wc_i[(wc_j < 0)|(wc_j >= hr_height)|(wc_i < 0)|(wc_i >= hr_width)] = -1

# map the landmask to the UTM grid
# bin_indices is the flattened wc coordinates mapped to hr grid bin indices
# landmask_values is the corresponding flattened wc values, 0 or 1
# majority is a function that returns 0 (water) for an empty array and else the majority of 0 or 1, 1 if equal
# hr_size is the size of the hr grid, i.e. the number of bins
# bin_values is the flattened array of values 0 water or 1 land in the hr grid
bin_indices = (wc_j * hr_width + wc_i).flatten()
del wc_i, wc_j
landmask_values = (landmask_wgs84.data).flatten()
# scipy is not the fastest, but it is convenient to do the resampling
bin_values, _edges, _indices = scipy.stats.binned_statistic(bin_indices, landmask_values, majority, hr_size, (0, hr_size))
scipy.stats.binned_statistics_dd()
del landmask_values, _edges, _indices
landmask = bin_values.reshape(hr_height, hr_width).astype(np.int8)
"""
watermask_da = xr.DataArray(landmask.reshape((1, hr_height, hr_width)), coords=dataset.coords, dims=dataset.dims, attrs=dataset.attrs)
watermask_da.rio.to_raster('worldcover-32UME.tif')
"""
del bin_values

# ----------------------------------------

print("mask inland water ...")

# add inland water in coastal or land areas
mask[((mask==LAND)|(mask==COASTAL))&(landmask==0)] = INLANDWATER
landmask = None

# ----------------------------------------

print("fix coastal 'inland water' ...")

# correct coastal "inland water" by extending ocean (64) into inland water (96) as "coastal inland water" (1) ...
depth = 44
TEMPORARY = 1
for i in range(depth):
    dilate(mask, OCEAN, four_neighbours, 1, INLANDWATER, TEMPORARY)
# ... and by reverting "costal inland water (1) to inland water (96) in river mouths
for i in range(depth):
    dilate(mask, INLANDWATER, four_neighbours, 1, TEMPORARY, INLANDWATER)
# replace non-reverted "coastal inland water" (1) by coastal land (160)
mask[mask==TEMPORARY] = COASTAL

# ----------------------------------------

print("mask transition zone from ocean to inland water ...")

# values 64 < v < 96 are the transition zone from ocean to inland water
# incremental values 7+distance from ocean up to 32 pixels
# 4 neighbours for the first two steps to avoid jumps over locks
depth = 32
for i in range(depth):
    dilate(mask, OCEAN+i, four_neighbours if i<2 else eight_neighbours, 1, INLANDWATER, OCEAN+i+1)

# ----------------------------------------

# dilate from transition zone into coastal and from inland water into coastal
# temporarily label shore of both
shore_width = 4
transition_depth = 32
SHORE_TRANSITION = 2
SHORE_INLANDWATER = 1
for i in range(shore_width):
    if i == 0:
        dilate_range(mask, OCEAN+1, INLANDWATER-1, eight_neighbours, 1, COASTAL, SHORE_TRANSITION)
    else:
        dilate(mask, SHORE_TRANSITION, eight_neighbours, 1, COASTAL, SHORE_TRANSITION)
    dilate(mask, INLANDWATER, eight_neighbours, 1, COASTAL, SHORE_INLANDWATER)
# count distance from ocean in shore of transition zone
for i in range(transition_depth):
    dilate(mask, OCEAN if i == 0 else COASTAL+i, eight_neighbours, 1, SHORE_TRANSITION, COASTAL+i+1)
# revert inland water shore and unreached transition zone shore
mask[(mask==SHORE_INLANDWATER)|(mask==SHORE_TRANSITION)] = COASTAL

# ----------------------------------------

print("write s2w mask ...")

s2w_mask_path = f"s2w-mask-{granule}.tif"
# write result to s2w mask
s2w = xr.DataArray(mask.reshape((1, hr_height, hr_width)), coords=dataset.coords, dims=dataset.dims, attrs=dataset.attrs)
s2w.rio.to_raster(s2w_mask_path, compress='LZW', tiled=True)

print(s2w_mask_path)

