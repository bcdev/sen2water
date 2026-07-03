# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

from typing import Tuple

# changes in 1.1:
# ...

import numpy as np
import dask.array as da
import xarray as xr
from eopf.computing import EOProcessingUnit, MappingDataType
from pyproj import Transformer


class GeoCoordinatesPU(EOProcessingUnit):

    def run(
        self,
        inputs: MappingDataType,
        transformer: Transformer = None,
        dtype: str = None,
        dims: Tuple[int, int] = None,
        chunks: Tuple[int, int] = None,
        **kwargs,
    ) -> MappingDataType:

        # check input

        if not "metric_coords" in inputs:
            raise KeyError("No 'metric_coords' in input of GeoCoordinatesPU.")
        metric_coords_dt: xr.DataTree = inputs["metric_coords"]
        if not isinstance(metric_coords_dt, xr.DataTree):
            raise TypeError("Input 'metric_coords' of GeoCoordinatesPU is not an xarray.DataTree.")

        result_data = da.map_blocks(
            self.geo_coordinates,
            metric_coords_dt["y"].data,
            metric_coords_dt["x"].data,
            transformer=transformer,
            new_axis=0,
            chunks=(2, *chunks),
            dtype=dtype,
            meta=np.array((), dtype=dtype),
            **kwargs,
        )

        result = xr.DataTree()
        result["latitude"] = xr.DataArray(
            result_data[0],
            dims=dims,
            attrs={"standard_name": "latitude",
                   "units": "degrees_north"},
        )
        result["longitude"] = xr.DataArray(
            result_data[1],
            dims=dims,
            attrs={"standard_name": "longitude",
                   "units": "degrees_east"},
        )
        return {"geographic_coords": result}

    def geo_coordinates(
        self,
        x,
        y,
        *,
        transformer: Transformer = None,
    ) -> np.ndarray:
        """
        Transforms UTM coordinates into geographic coordinates

        x: np.ndarray
        """

        lat, lon = transformer.transform(x, y)
        return np.stack((lat, lon))

    func = geo_coordinates
