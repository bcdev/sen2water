#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Generate Sen2Water mask"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "proprietary, TBD"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

print("loading software ... ")

import rioxarray as rio
import xarray as xr
import numpy as np
import dask.array as da
import dask
from multiprocessing.pool import ThreadPool
import scipy
from pyproj import CRS, Transformer
from pyproj.enums import TransformDirection
import math
import sys

# define some constants and functions

INVALID = 0
OCEAN = 64
INLANDWATER = 96
BOTTOMREFLECTION = 128
COASTAL = 160
LAND = 192

four_neighbours=np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=np.uint8)
eight_neighbours=np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=np.uint8)

def majority(a):
    return (np.sum(a) + len(a) // 2) // len(a) if len(a) > 0 else 0

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

def calc_extent(hr_y, hr_x, trafo=None, global_rows=None, global_cols=None, hr_step=None, block_id=None):
    half_a_pixel = hr_step / 2
    # extract image border (top, bottom, left, right) for transformation
    hr_y_border = np.stack((hr_y[0] + half_a_pixel,
                            hr_y[-1] - half_a_pixel,
                            hr_y[:,0],
                            hr_y[:,-1]))
    hr_x_border = np.stack((hr_x[0],
                            hr_x[-1],
                            hr_x[:,0] - half_a_pixel,
                            hr_x[:,-1] + half_a_pixel))
    # move four pixels of border into image corners
    hr_x_border[0,0] -= half_a_pixel
    hr_x_border[0,-1] += half_a_pixel
    hr_x_border[1,0] -= half_a_pixel
    hr_x_border[1,-1] += half_a_pixel
    # transform border from UTM to WGS84 geographic coordinates
    hr_lat,hr_lon = trafo.transform(hr_x_border.flatten(), hr_y_border.flatten(), direction=TransformDirection.INVERSE)
    # determine minmax lat of border
    hr_lat_min = np.min(hr_lat)
    hr_lat_max = np.max(hr_lat)
    # turn lon by first lon as reference to determine minmax lon of border
    hr_lon_ref = hr_lon[0]
    hr_lon_min = (np.min((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
    hr_lon_max = (np.max((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
    # determine minmax pixel coordinates of worldcover subset
    wc_j_min = int(math.floor((90.0 - hr_lat_max) / 180.0 * global_rows))
    wc_j_max = int(math.ceil((90.0 - hr_lat_min) / 180.0 * global_rows))
    wc_i_min = int(math.floor((hr_lon_min + 180.0) / 360.0 * global_cols))
    wc_i_max = int(math.ceil((hr_lon_max + 180.0) / 360.0 * global_cols))
    i_step = int(1.0 / math.cos((hr_lat_min + hr_lat_max) / 2 * math.pi / 180.0))
    #print(f"extent {block_id=} {wc_j_min=} {wc_j_max=} {wc_i_min=} {wc_i_max=}")
    return np.array([[[wc_j_min]], [[wc_j_max]], [[wc_i_min]], [[wc_i_max]], [[i_step]]])

def globalmask_subset(global_mask, j_min, j_max, i_min, i_max, i_step):
    if i_min <= i_max:
        mask_subset = global_mask[j_min:j_max, i_min:i_max:i_step]
    else:
        mask_subset = xr.concat((global_mask[j_min:j_max, i_min:], global_mask[j_min:j_max, :i_max]), dim='x')
        mask_subset = mask_subset[:,::i_step]
    return mask_subset

def globalmask_bincounts(bin_indices, block_size):
    return np.bincount(bin_indices, minlength=block_size)

def globalmask_binning(mask_values, bin_indices, all_bin_counts, block_height, block_width):
    # map the landmask to the UTM grid
    # bin_indices is the flattened wc coordinates mapped to hr grid bin indices
    # mask_values is the corresponding flattened wc values, 0 or 1
    block_size = block_height * block_width
    land_bin_counts = np.bincount(bin_indices[mask_values != 0], minlength=block_size)
    bin_values = 2 * land_bin_counts >= all_bin_counts
    del land_bin_counts
    return bin_values.reshape((block_height, block_width)).astype(np.uint8)

def target_pixel_indices(lat_vector, lon_vector, trafo, block_top, block_left, hr_step):
    # transform land mask coordinates from geographic WGS84 to metric UTM
    lat = np.tile(lat_vector.values, (len(lon_vector), 1)).T
    lon = np.tile(lon_vector.values, (len(lat_vector), 1))
    x, y = trafo.transform(lat, lon, direction=TransformDirection.FORWARD)
    del lat, lon
    # determine block pixel coordinates in HR-OC UTM grid for landmask subset
    j = ((block_top - y + hr_step / 2) / hr_step).astype(np.int32)
    i = ((x - block_left + hr_step / 2) / hr_step).astype(np.int32)
    return j, i

def calc_combined_mask(target_block, extents,
                       trafo=None, worldcover_path=None, globalislands_path=None, continentalshoreline_path=None,
                       hr_step=0, hr_left=0, hr_top=0,
                       block_info=None, block_id=None):
    wc_j_min = extents[0,0,0]
    wc_j_max = extents[1,0,0]
    wc_i_min = extents[2,0,0]
    wc_i_max = extents[3,0,0]
    i_step = extents[4,0,0]
    # block metric coordinates and extent
    block_top = hr_top - block_info[0]['array-location'][0][0] * hr_step
    block_left = hr_left + block_info[0]['array-location'][1][0] * hr_step
    block_height, block_width = target_block.shape
    block_size = block_height * block_width
    # read subset, determine mapping from subset to target pixels, determine mask of subset pixels inside target, perform binning
    print("determine mapping from geographic to UTM grid ...")
    with xr.open_dataset(worldcover_path, engine='zarr', mask_and_scale=False) as wc:
        wc_subset = globalmask_subset(wc.landmask, wc_j_min, wc_j_max, wc_i_min, wc_i_max, i_step)
        wc_j, wc_i = target_pixel_indices(wc_subset.y, wc_subset.x, trafo, block_top, block_left, hr_step)
        # exclude pixels outside of HR-OC mask granule
        is_inside_granule = ((wc_j >= 0) & (wc_j < block_height) & (wc_i >= 0) & (wc_i < block_width)).astype(bool)
        bin_indices = (wc_j * block_width + wc_i)[is_inside_granule]
        all_bin_counts = globalmask_bincounts(bin_indices, block_size)
        #print(f"{is_inside_granule.shape=} {len(bin_indices)=}")
        del wc_j, wc_i
        print("subset and reproject worldcover land mask ...")
        wc_landmask = globalmask_binning(wc_subset.data[is_inside_granule], bin_indices, all_bin_counts, block_height, block_width)
        del wc_subset
    # mark as land what is land in worldcover
    mask = wc_landmask
    del wc_landmask
    mask[(mask != 0)] = LAND
    # read subset of globalislands, perform binning
    print("subset and reproject globalislands land mask ...")
    with xr.open_dataset(globalislands_path, engine='zarr', mask_and_scale=False) as gi:
        gi_subset = globalmask_subset(gi.inlandmask, wc_j_min, wc_j_max, wc_i_min, wc_i_max, i_step)
        gi_landmask = globalmask_binning(gi_subset.data[is_inside_granule], bin_indices, all_bin_counts, block_height, block_width)
        del gi_subset
    # mark as inland water what is not land in worldcover and not ocean in globalislands
    mask[((mask == 0) & (gi_landmask!=0))] = INLANDWATER
    del gi_landmask
    # read subset of continentalshoreline, perform binning
    print("subset and reproject continentalshoreline land mask ...")
    with xr.open_dataset(continentalshoreline_path, engine='zarr', mask_and_scale=False) as cs:
        cs_subset = globalmask_subset(cs.continentalmask, wc_j_min, wc_j_max, wc_i_min, wc_i_max, i_step)
        cs_landmask = globalmask_binning(cs_subset.data[is_inside_granule], bin_indices, all_bin_counts, block_height, block_width)
        del cs_subset
    del bin_indices, all_bin_counts, is_inside_granule
    # mark as ocean what is not land and not inland water and is ocean in continental shoreline
    mask[((mask == 0) & (cs_landmask==0))] = OCEAN
    del cs_landmask
    # mark as river water the remaining rivers connected to ocean not marked as land by globalislands
    RIVERWATER=3
    mask[(mask==0)] = RIVERWATER
    #
    print(f"mask {block_id} done")
    return mask

def buffer_ocean_from_continentalshoreline_to_worldcover(mask):
    RIVERWATER=3
    depth = 24  # 32
    # buffer ocean into rivers (to remove non-river pixels) ...
    for i in range(depth):
        dilate(mask, OCEAN, four_neighbours, 1, RIVERWATER, OCEAN)
    # ... and revert in river mouths
    for i in range(depth):
        dilate(mask, RIVERWATER, four_neighbours, 1, OCEAN, RIVERWATER)
    mask[mask==RIVERWATER] = INLANDWATER

def buffer_ocean_from_globalislands_to_worldcover(mask):
    # Are there lakes close to shoreline that are partially marked as ocean?
    # buffer ocean into inlandwater
    depth = 48  # 24
    for i in range(depth):
        dilate(mask, OCEAN, four_neighbours, 1, INLANDWATER, OCEAN)
    # buffer inlandwater back into ocean in estuaries
    for i in range(depth):
        dilate(mask, INLANDWATER, four_neighbours, 1, OCEAN, INLANDWATER)

def buffer_ocean_into_land_as_coastal(mask):
    depth = 28
    for i in range(depth):
        dilate(mask, OCEAN, eight_neighbours, 1, LAND, COASTAL)

def label_coastal_transition_zone(mask):
    # dilate from transition zone into coastal and from inland water into coastal
    # temporarily label shore of both
    # values 64 < v < 96 are the transition zone from ocean to inland water
    # incremental values 7+distance from ocean up to 32 pixels
    # 4 neighbours for the first two steps to avoid jumps over locks
    depth = 32
    for i in range(depth):
        dilate(mask, OCEAN+i, four_neighbours if i<2 else eight_neighbours, 1, INLANDWATER, OCEAN+i+1)
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


if __name__ == '__main__':
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
    countslist_path = "/data/sen2water/globalmaskcounts.list"
    zeromask_path = "zero-masks/01CDJ-zeros.tif"
    s2w_mask_path = "global-masks/s2w-globalmask-01CDJ.tif"
    zeromask_path = "zero-masks/32UME-zeros.tif"
    s2w_mask_path = "global-masks/s2w-globalmask-32UME.tif"
    """

    worldcover_path = sys.argv[1]
    globalislands_path = sys.argv[2]
    continentalshoreline_path = sys.argv[3]
    countslist_path = sys.argv[4]
    zeromask_path = sys.argv[5]
    s2w_mask_path = sys.argv[6]

    chunksize = 610
    
    dask.config.set(scheduler="threads")
    dask.config.set(pool=ThreadPool(9))
    #dask.config.set(scheduler="synchronous")
    #from dask.distributed import Client, LocalCluster
    #import time
    #cluster = LocalCluster(processes=True, n_workers=9, threads_per_worker=1)
    #client = Client(cluster)
    #print(client.dashboard_link)
    #time.sleep(10)

    #import cProfile
    #profile = cProfile.Profile()
    #profile.enable()

    print("determine extents of worldcover grid for granule in UTM grid ...")
    
    dataset = rio.open_rasterio(zeromask_path, chunks={"y": chunksize, "x": chunksize})
    hr_width = len(dataset.x)
    hr_height = len(dataset.y)
    hr_size = hr_height * hr_width
    hr_step = 60.0
    hr_left = dataset.x.data[0]
    hr_top = dataset.y.data[0]
    
    #wc = xr.open_dataset(worldcover_path, engine='zarr', mask_and_scale=False)
    #gi = xr.open_dataset(globalislands_path, engine='zarr', mask_and_scale=False)
    #cs = xr.open_dataset(continentalshoreline_path, engine='zarr', mask_and_scale=False)
    
    with xr.open_dataset(worldcover_path, engine='zarr', mask_and_scale=False) as wc:
        # create geo-trafo WGS84 -> UTM
        hr_crs = CRS.from_cf(dataset.spatial_ref.attrs)
        wc_crs = CRS.from_cf(wc.spatial_ref.attrs)
        wgs_to_utm = Transformer.from_crs(wc_crs, hr_crs)
        
        hr_y = da.from_array(np.tile(dataset.y.values, (hr_width, 1)).T, chunks=dataset.data[0].chunks)
        hr_x = da.from_array(np.tile(dataset.x.values, (hr_height, 1)), chunks=dataset.data[0].chunks)

    extents = da.map_blocks(calc_extent,
                            hr_y, hr_x,
                            trafo=wgs_to_utm, global_rows=wc.y.shape[0], global_cols=wc.x.shape[0], hr_step=hr_step,
                            new_axis=0, chunks=(4,1,1), dtype=np.int32)
    
    combined_mask = da.map_blocks(calc_combined_mask,
                                  dataset.data[0],
                                  extents,
                                  trafo=wgs_to_utm,
                                  worldcover_path=worldcover_path, globalislands_path=globalislands_path, continentalshoreline_path=continentalshoreline_path,
                                  hr_step=hr_step, hr_left=hr_left, hr_top=hr_top,
                                  drop_axis=0, chunks=dataset.chunks[1:], dtype=np.uint8)

    #print("starting computation")
    
    mask = dask.compute(combined_mask)[0]
    
    print("buffering ocean into rivers and back to preserve ocean included in continental shoreline ...")
    buffer_ocean_from_continentalshoreline_to_worldcover(mask)
    print("buffering ocean into inlandwater and back to extend ocean up to worldcover coastline ...")
    buffer_ocean_from_globalislands_to_worldcover(mask)
    print("buffering ocean into land to mark coastal land ...")
    buffer_ocean_into_land_as_coastal(mask)
    print("mask transition zone from ocean to inland water ...")
    label_coastal_transition_zone(mask)

    print("counting pixel classes ...")

    land_count = np.count_nonzero(mask>=128)
    inlandwater_count = np.count_nonzero((mask<128) & (mask > 64))
    ocean_count = np.count_nonzero(mask<=64)
    with open(countslist_path, "a") as f:
        f.write(f"{s2w_mask_path}\t{ocean_count}\t{inlandwater_count}\t{land_count}\n")

    print("write s2w mask ...")
    s2w = xr.DataArray(mask.reshape((1, hr_height, hr_width)), coords=dataset.coords, dims=dataset.dims, attrs=dataset.attrs)
    s2w.rio.to_raster(s2w_mask_path, compress='LZW', tiled=True)

    #profile.disable()
    #import io
    #import pstats
    #buffer = io.StringIO()
    #stats = pstats.Stats(profile, stream=buffer)
    #stats.sort_stats("tottime").print_stats()
    #with open("prof.out", "w") as f:
    #    f.write(buffer.getvalue())

    #client.close()
    #cluster.close()

    print(s2w_mask_path)
