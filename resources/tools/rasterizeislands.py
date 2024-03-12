# -*- coding: utf-8 -*-

"""
Tool to ...
"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2024, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

from glob import glob
from pmonitor import PMonitor

# islands-by-tile/shoreline-h059v10.json
input_root = 'islands-per-tile'
output_root = 'islands-raster'

inputs = glob(f"{input_root}//*json')
pm = PMonitor(inputs,
              request='rasterizeislands',
              hosts=[('localhost', 8)],
              types=[('gdal_rasterize', 8)],
              logdir='log',
              simulation=False)
for input in inputs:
    row = int(input[31:33])
    col = int(input[27:30])
    lat1 = 90 - row * 3
    lat2 = lat1 - 3
    lon1 = col * 3 - 180
    lon2 = lon1 + 3
    # gdal_rasterize -l GlobalIslandsv2_WGS84_fixed -burn 1 -te -3 57 0 60 -ts 36000 36000 \
    #     -ot uint8 -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 \
    #     -co NUM_THREADS 4 islands-by-tile/shoreline-h059v10.json globalislands-h059v10.tif
    pm.execute('gdal_rasterize',
               [ input ],
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

