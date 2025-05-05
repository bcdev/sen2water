#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Converts GSHHS into the WorldCover raster"""

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

# islands-by-tile/shoreline-h059v10.json
coastal_tiles = "WorldCover_coastal_tiles.geojsonl.json"
land_tiles = "WorldCover_land_tiles.geojsonl.json"
inputs = [ "GSHHS_L1_gt10000sqkm_buffered_minus033_edit.json" ]
output_root = 'continentalshoreline-byte-raster'

tiles = []
try:
    for list in [coastal_tiles, land_tiles]:
        with open(list) as f:
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

pm = PMonitor(inputs,
              request='continentalshoreline',
              hosts=[('localhost', 8)],
              types=[('gdal_rasterize', 8)],
              logdir='log',
              simulation=False)
for tile in tiles:
    row, col = tile
    lat1 = 90 - row * 3
    lat2 = lat1 - 3
    lon1 = col * 3 - 180
    lon2 = lon1 + 3
    # gdal_rasterize -l GSHHS_L1_gt10000sqkm_buffered_minus033_edit -burn 1 -te -3 57 0 60 \
    #     -ts 36000 36000 -ot Byte -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 \
    #     -co NUM_THREADS 4 GSHHS_L1_gt10000sqkm_buffered_minus033_edit.json continentalshoreline-h059v10.tif
    pm.execute('gdal_rasterize',
               inputs,
               [ f"{output_root}/continentalshoreline-h{col:03d}v{row:02d}.tif" ],
               [ "-l", "GSHHS_L1_gt10000sqkm_buffered_minus033_edit",
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
