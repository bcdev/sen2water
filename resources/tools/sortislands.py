# -*- coding: utf-8 -*-

"""
Tool to sort islands represented by polygons in GeoJSON into coastal 3 degree worldcover tiles.

The tool makes use of the fact that coastal tiles are not included completely in a continent.
Only edges of polygons overlapping with the tile need to be considered.
The method identifies all tiles of the bounding box of each edge of the polygon of the island.
Multipolygons are handled by considering each polygon included.
An island may contribute to more than a single tile.
"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "proprietary, TBD"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import json

# ogr2ogr -f GeoJSON globalIslandsv2_WGS84_fixed.json GlobalIslandsv2_WGS84_fixed.shp
coastal_tiles = "WorldCover_coastal_tiles.geojsonl.json"
global_islands = "globalIslandsv2_WGS84_fixed.json"
island_rasters_root = "islands-by-tile"

# collect coastal tile index

resolution = 3
tiles = {
#    (10,59): None,
#    (9,59): None,
}
try:
  with open(coastal_tiles) as f:
    while True:
        line = f.readline()
        if line is None or len(line) == 0: break
        record = json.loads(line[:-1])
        lon, lat = record['geometry']['coordinates'][0][0][0]
        tile = (90 - int(lat)) // 3, (180 + int(lon)) // 3
        tiles[tile] = None
except Exception as e:
    print(f"{e=}")
    print(f"{line=}")
    raise e

# functions to sort edges of polygons by tile index

def collect_row(j, i1, i2, line_tiles):
    if i1 <= i2:
        if i2 - i1 <= 60:
            # common eastward edge in direction of positive i
            for i in range(i1, i2+1):
                line_tiles.add((j, i))
        else:
            # antimeridian-crossing westward edge
            for i in range(i2, i1 + 120 + 1):
                line_tiles.add((j, i % 120))
    else:
        if i1 - i2 <= 60:
            # common westward edge in direction of negateve i
            for i in range(i2, i1+1):
                line_tiles.add((j, i))
        else:
            # antimeridian-crossing eastward edge
            for i in range(i1, i2 + 120 + 1):
                line_tiles.add((j, i % 120))

def collect_line_tiles(polygon, resolution, line_tiles):
    lon1, lat1 = None, None
    for lon, lat in polygon:
        if lon1 is not None:
            # edge from lat1/lon1 to lat/lon
            # transform into worldcover tile grid
            j1 = int((90 - lat1) / 3)
            i1 = int((180 + lon1) / 3)
            j2 = int((90 - lat) / 3)
            i2 = int((180 + lon) / 3)
            if j1 <= j2:
                # southward edge in direction of positive j
                for j in range(j1, j2+1):
                    collect_row(j, i1, i2, line_tiles)
            else:
                # northward edge in direction of negative j
                for j in range(j2, j1+1):
                    collect_row(j, i1, i2, line_tiles)
        lon1, lat1 = lon, lat

# read hand process global islands island polygon by island polygon

try:
  with open(global_islands) as f:
    header = [ f.readline() for i in range(5) ]
    c = 5
    while True:
        line = f.readline()
        if line is None: break
        if line == ']\n': break
        # collect all tile indices touched by edges of the polygon
        line_tiles = set()
        island = json.loads(line[:-2 if line[-2] == ',' else -1])
        id = island["properties"]["OBJECTID_1"]
        name = island["properties"]["NAME_LOCAL"]
        polygon_type = island['geometry']['type']
        coordinates = island['geometry']['coordinates']
        if polygon_type == 'Polygon':
            for polygon in coordinates:
                collect_line_tiles(polygon, resolution, line_tiles)
        elif polygon_type == 'MultiPolygon':
            for multi in coordinates:
                for polygon in multi:
                    collect_line_tiles(polygon, resolution, line_tiles)
        else:
            raise ValueError("unknown polygon type {polygon_type}")
        # distribute polygon to all coastal tiles touched by the polygon
        c2 = 0
        for tile in line_tiles:
            if tile in tiles:
                # lazily create polygon file for a tile
                if tiles[tile] is None:
                    (j, i) = tile
                    tiles[tile] = open(f"{island_rasters_root}/shoreline-h{i:03d}v{j:02d}.json", "w")
                    for l in header:
                        tiles[tile].write(l)
                tiles[tile].write(line)
                print(f"adding line {c} {name} to {tile}")
                c2 += 1
        c += 1
    # close coastal tiles files
    for tile in tiles:
        if tiles[tile] is not None:
            tiles[tile].write("]}\n")
            tiles[tile].close()
except Exception as e:
    print(e)
    print(f"{c=}")
    print(f"{id=}")
    print(f"{polygon=}")
    
print("done")

