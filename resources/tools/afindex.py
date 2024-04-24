#!/usr/bin/env python

print("loading software ... ")

import glob
import math
import rioxarray as rio
import numpy as np
from pyproj import CRS, Transformer
from pyproj.enums import TransformDirection
import sys

if len(sys.argv) != 5:
    print('Usage:')
    print('  afindex.py <affile> <countslist> <zerogranulepath> <outputpath>')
    print('e.g.')
    print("  afindex.py viirs-firms-2023.csv afcounts.list zero-masks af-index")
    sys.exit(1)

"""
af_path = "/data/fire/viirs-firms-2023.csv"
countslist_path = "afcounts.list"
zeromask_path = "zero-masks"
index_path = "af-index"
template = "zero-masks/19FEA-zeros.tif"
"""

af_path = sys.argv[1]
countslist_path = sys.argv[2]
zeromask_path = sys.argv[3]
index_path = sys.argv[4]

print("reading active fires ...")

height = 180
width = 360
step = 1.0
geoindex = dict()

with open(af_path) as f:
    line = f.readline()
    count = 1
    while True:
        line = f.readline()
        if not line: break
        count += 1
        try:
            fields = line.split(',')
            lat = float(fields[0])
            lon = float(fields[1])
            date = fields[5]  # TODO which field is date?
            i = int((90.0 - lat) / step) * width + int((lon + 180.0) / step)
            if i in geoindex:
                geoindex[i].append((lat, lon, date))
            else:
                geoindex[i] = [(lat, lon, date)]
        except Exception as e:
            print(f"error {e} in line {count}: {line}")
            raise e

for i in geoindex:
    lat = np.array([entry[0] for entry in geoindex[i]])
    lon = np.array([entry[1] for entry in geoindex[i]])
    date = np.array([entry[2] for entry in geoindex[i]])
    geoindex[i] = (lat, lon, date)

print(f"{len(geoindex)=}")
print("reading template granules ...")

"""
if True:
"""
for template in glob.glob(f"{zeromask_path}/*tif"):
    granule = template[template.rfind('/')+1:template.rfind('/')+6]

    with rio.open_rasterio(template) as hr:
        hr_width = len(hr.x)
        hr_height = len(hr.y)
        hr_size = hr_height * hr_width
        hr_step = 60.0
        hr_left = hr.x.data[0]
        hr_top = hr.y.data[0]

        # create geo-trafo WGS84 -> UTM
        hr_crs = CRS.from_cf(hr.spatial_ref.attrs)
        af_crs = CRS.from_epsg(4326)
        wgs_to_utm = Transformer.from_crs(af_crs, hr_crs)

        # map granule to geographic grid
        hr_y = hr.y.values
        hr_x = hr.x.values
        hr_y_m = np.tile(hr_y, (hr_x.shape[0],1)).T
        hr_x_m = np.tile(hr_x, (hr_y.shape[0],1))
        # transform border from UTM to WGS84 geographic coordinates
        hr_lat,hr_lon = wgs_to_utm.transform(hr_x_m.flatten(), hr_y_m.flatten(), direction=TransformDirection.INVERSE)

        # determine bounding box in geographic grid
        hr_lat_min = np.min(hr_lat)
        hr_lat_max = np.max(hr_lat)
        # turn lon by first lon as reference to determine minmax lon of border
        hr_lon_ref = hr_lon[0]
        hr_lon_min = (np.min((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
        hr_lon_max = (np.max((hr_lon - hr_lon_ref + 180.0) % 360.0 - 180.0) + hr_lon_ref + 180.0) % 360.0 - 180.0
        # determine minmax pixel coordinates of index subset
        af_j_min = int(math.floor((90.0 - hr_lat_max) / 180.0 * height))
        af_j_max = int(math.ceil((90.0 - hr_lat_min) / 180.0 * height))
        af_i_min = int(math.floor((hr_lon_min + 180.0) / 360.0 * width))
        af_i_max = int(math.ceil((hr_lon_max + 180.0) / 360.0 * width))
        print(f"extent {af_j_min=} {af_i_min=} {af_j_max-af_j_min=} {af_i_max-af_i_min=}")
        if hr_lon_min > hr_lon_max:
            print('crossing antimeridian')

        dates = []
        for j in range(af_j_min, af_j_max + 1):
            if af_i_min <= af_i_max:
                for i in range(af_i_min, af_i_max + 1):
                    ji = j * width + i
                    if ji in geoindex:
                        lat, lon, date = geoindex[ji]
                        x, y = wgs_to_utm.transform(lat, lon, direction=TransformDirection.FORWARD)
                        print(date[(x >= hr_x[0]) & (x <= hr_x[-1]) & (y >= hr_y[-1]) & (y <= hr_y[0])])
                        if len(date) > 0:
                            print(date)
                        dates.extend(list(date[(x >= hr_x[0]) & (x <= hr_x[-1]) & (y >= hr_y[-1]) & (y <= hr_y[0])]))
            else:
                for i in range(af_i_min, width + af_i_max + 1):
                    ji = j * width + (i % width)
                    if ji in geoindex:
                        lat, lon, date = geoindex[ji]
                        x, y = wgs_to_utm.transform(lat, lon, direction=TransformDirection.FORWARD)
                        dates.extend(list(date[(x >= hr_x[0]) & (x <= hr_x[-1]) & (y >= hr_y[-1]) & (y <= hr_y[0])]))

    print(f"granule {granule} has {len(dates)} af")

    with open(countslist_path, "a") as f:
        f.write(f"{granule}\t{len(dates)}\n")
    if len(dates) > 0:
        with open(f"{index_path}/{granule}.dates", "w") as f:
            for d in sorted(list(dates)):
                f.write(f"{d}\n")
        print(f"{index_path}/{granule}.dates")

print("done")
