from typing import Any, Optional, Literal, Tuple
from collections.abc import Mapping

import xarray as xr
import dask
import dask.array as da
from pyproj import Transformer, CRS

from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from sen2water.msiresampling.geocoordinates import GeoCoordinates


class ResamplingPU(EOProcessingUnit):

    PROCESSOR_NAME = "MsiResampling"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    INPUT_IDENTIFIER: str = "l1c"

    def __init__(self, identifier: Any = ""):
        super().__init__(identifier)
        self._update_history = True

    def run(
        self,
        inputs: MappingDataType,
        adfs: Optional[MappingAuxiliary] = None,
        *,
        resolution: int = 60,
        chunksize_in_meters = 36600,
        downsampling: Literal['detectormean', 'first', 'min', 'max', 'mean', 'median'] = "detectormean",
        flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ] = "flagand",
        upsampling: Literal["nearest", "bilinear", "bicubic"] = "nearest",
        with_detfoo_filter: bool = False,
        merge_flags: bool = False,

        **kwargs
    ) -> MappingDataType:
        if len(inputs) < 1:
            raise MyError("Missing mandatory input product 'l1c'")

        l1c = inputs[self.INPUT_IDENTIFIER]

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {"y": l1c.sizes[f"y{resolution}"], "x": l1c.sizes[f"x{resolution}"]}
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)
        # TODO necessary?
        resampled = xr.DataTree(name="resampled")
        coordinate_bands = {}
        target_bands = {}
        #resampled["measurements/reflectance/resampled/x"]


        # copy 1-d metric coordinate bands

        coordinate_bands["y"] = xr.DataArray(
            l1c[f"measurements/reflectance/r{resolution}m/y"].data,
            dims=['y'],
            attrs={"long_name": "y coordinate of projection",
                   "standard_name": "projection_y_coordinate",
                   "units": "m"}
        )
        coordinate_bands["x"] = xr.DataArray(
            l1c[f"measurements/reflectance/r{resolution}m/x"].data,
            dims=['x'],
            attrs={"long_name": "x coordinate of projection",
                   "standard_name": "projection_x_coordinate",
                   "units": "m"}
        )
        yy, xx = self.construct_coordinate_arrays(
            l1c[f"measurements/reflectance/r{resolution}m/y"],
            l1c[f"measurements/reflectance/r{resolution}m/x"],
            chunks=chunks,
        )

        transformer = Transformer.from_crs(CRS.from_epsg(l1c.attrs["stac_discovery"]["properties"]["proj:code"]), CRS("EPSG:4326"))
        lat_lon = GeoCoordinates().apply(
            xx,
            yy,
            transformer=transformer,
            new_axis=0,
            dtype=np.float64,
            chunks=(2, *chunks),
        )
        lat_data = lat_lon[0]
        lon_data = lat_lon[1]
        del x, y, xx, yy, transformer, lat_lon

    @dask.delayed
    def construct_coordinate_arrays(self, y, x, chunks: Tuple[int, int]):
        y = da.from_array(y, chunks=chunks[0])
        x = da.from_array(x, chunks=chunks[1])
        xx = da.broadcast_to(x, (y.shape[0], x.shape[0]), chunks=chunks)
        yy = da.transpose(da.broadcast_to(y, (y.shape[0], x.shape[0]), chunks=chunks[::-1]))
        return yy, xx



