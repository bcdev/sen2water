#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Control script to rasterize GlobalIslands"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "proprietary, TBD"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import json
from pmonitor import PMonitor

inland_tiles = "WorldCover_land_tiles.geojsonl.json"
global_geojson = "global.json"
output_root = "land-byte-rasters"

pm = PMonitor([global_geojson],
              request='globalinland',
              hosts=[('localhost', 8)],
              types=[('gdal_rasterize', 8)],
              logdir='log',
              simulation=False)

tiles = []
try:
  with open(inland_tiles) as f:
    while True:
        line = f.readline()
        if line is None or len(line) == 0: break
        record = json.loads(line[:-1])
        lon, lat = record['geometry']['coordinates'][0][0][0]
        tile = (90 - int(lat)) // 3, (180 + int(lon)) // 3
        tiles.append(tile)
except Exception as e:
    print(f"{e=}")
    print(f"{line=}")
    raise e

for tile in tiles:
    row, col = tile
    lat1 = 90 - row * 3
    lat2 = lat1 - 3
    lon1 = col * 3 - 180
    lon2 = lon1 + 3
    # gdal_rasterize -l GlobalIslandsv2_WGS84_fixed -burn 1 -te -3 57 0 60 -ts 36000 36000 \
    #     -ot uint8 -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 \
    #     -co NUM_THREADS 4 global.json globalislands-h059v10.tif
    pm.execute('gdal_rasterize',
               [ global_geojson ],
               [ f"{output_root}/globalislands-h{col:03d}v{row:02d}.tif" ],
               [ "-l", "GlobalIslandsv2_WGS84_fixed",
                 "-burn", "1",
                 "-te", str(lon1), str(lat2), str(lon2), str(lat1),
                 "-ts", "36000", "36000",
                 "-ot", "Byte",
                 "-co", "COMPRESS=LZW",
                 "-co", "TILED=YES",
                 "-co", "BLOCKXSIZE=1024",
                 "-co", "BLOCKYSIZE=1024",
                 "-co", "NUM_THREADS=4" ])

pm.wait_for_completion()
