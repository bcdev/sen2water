from typing import Any, Optional, Literal, Tuple, Dict, List

import xarray as xr
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from sen2water.msiidepix.constants import IdepixMsiConstants as ic
from sen2water.msiidepix.idepixclassificationpu_pu import IdepixClassificationPUPU
from sen2water.msiidepix.idepixclassificationpu_algo import IdepixClassificationAlgoPU


class IdepixMainPU(EOProcessingUnit):
    PROCESSOR_NAME = "MsiIdepix"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    bands = ic.bands
    resolutions = ic.resolutions

    def __init__(self, identifier: Any = ""):
        super().__init__(identifier)
        self._update_history = True

    def run(
            self,
            inputs: MappingDataType,
            adfs: Optional[MappingAuxiliary] = None,
            mode: Optional[str] = None,
            *,
            resolution: int = 60,
            chunksize_in_meters=36600,
            ancillary: Optional[List[str]] = None,
            downsampling: Literal['detectormean', 'first', 'min', 'max', 'mean', 'median'] = "detectormean",
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ] = "flagand",
            upsampling: Literal["nearest", "bilinear", "bicubic"] = "nearest",
            with_detfoo_filter: bool = False,
            **kwargs
    ) -> MappingDataType:

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c: xr.DataTree = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        # determine parameters

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {
            "y": l1c[f"measurements/reflectance/resampled"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/resampled"].coords.sizes["x"]
        }
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)

        # create result accumulator

        idepix = xr.DataTree(name="idepix")

        # add geo-coding with 1-d metric coordinates, CRS, 2-d geographic coordinates

        # geographic_coords_dt = self._add_geocoding(
        #     l1c,
        #     idepix,
        #     ancillary,
        #     resolution,
        #     dims,
        #     chunks,
        # )

        # compute idepix...

        self._compute_idepix(
            l1c, idepix,
            flagdownsampling,
            resolution, chunksize_in_meters, dims, chunks
        )

        # copy coarse grained meteo data

        # idepix["conditions/meteorology/ecmwf"] = l1c["conditions/meteorology/ecmwf"].copy(inherit=False, deep=False)
        # idepix["conditions/meteorology/cams"] = l1c["conditions/meteorology/cams"].copy(inherit=False, deep=False)
        #
        # properties = l1c.attrs["stac_discovery"]["properties"]
        # idepix.attrs = {
        #     "other_metadata": {
        #         "product_type": "S2_MSI_Level-1C",
        #         "platform": properties["platform"],
        #         "start_date": properties["start_datetime"],
        #         "stop_date": properties["end_datetime"],
        #     }
        # }

        return {
            "idepix": idepix
        }

    # def _add_geocoding(
    #         self,
    #         l1c: xr.DataTree,
    #         idepix: xr.DataTree,
    #         ancillary: List[str] | None,
    #         resolution: int,
    #         dims: dict[str, int],
    #         chunks: tuple[int, int],
    # ) -> xr.DataTree:
    #
    #     epsg_code = l1c.attrs["stac_discovery"]["properties"]["proj:code"]
    #     # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
    #     # GT(1) w-e pixel resolution / pixel width.
    #     # GT(2) row rotation (typically zero).
    #     # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
    #     # GT(4) column rotation (typically zero).
    #     # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
    #     proj_transform = l1c.attrs["stac_discovery"]["properties"]["proj:transform"]
    #
    #     geo_transform = " ".join(str(proj_transform[i]) for i in [2, 0, 1, 5, 3, 4])
    #     crs = CRS.from_authority(*epsg_code.split(":"))
    #     transformer = Transformer.from_crs(crs, CRS("EPSG:4326"))
    #
    #     # The metric coordinates are numpy arrays, they are written as a single chunk.
    #     y_1d = l1c[f"measurements/reflectance/r{resolution}m/y"].data
    #     x_1d = l1c[f"measurements/reflectance/r{resolution}m/x"].data
    #     # Here we construct dask arrays to make the 2-d lat and lon chunked.
    #     y_data = da.from_array(y_1d, chunks=chunks[1])
    #     x_data = da.from_array(x_1d, chunks=chunks[0])
    #     y_grid = da.transpose(da.broadcast_to(y_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks[::-1]))
    #     x_grid = da.broadcast_to(x_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks)
    #
    #     metric_coords_dt = xr.DataTree()
    #     metric_coords_dt["y"] = xr.DataArray(y_grid, dims=dims, )
    #     metric_coords_dt["x"] = xr.DataArray(x_grid, dims=dims, )
    #
    #     geographic_coords_dt = GeoCoordinatesPU().run(
    #         {"metric_coords": metric_coords_dt},
    #         transformer=transformer,
    #         dtype=np.float64,
    #         dims=dims,
    #         chunks=chunks,
    #     )
    #
    #     y_da = xr.DataArray(
    #         y_1d,
    #         dims=["y"],
    #         attrs={
    #             "long_name": "y coordinate of projection",
    #             "standard_name": "projection_y_coordinate",
    #             "units": "m",
    #         },
    #     )
    #     x_da = xr.DataArray(
    #         x_1d,
    #         dims=["x"],
    #         attrs={
    #             "long_name": "x coordinate of projection",
    #             "standard_name": "projection_x_coordinate",
    #             "units": "m",
    #         },
    #     )
    #     crs_da = xr.DataArray(
    #         0,
    #         attrs={
    #             "crs_wkt": crs.to_wkt(),
    #             "GeoTransform": geo_transform
    #         }
    #     )
    #
    #     idepix["conditions/geometry/y"] = y_da
    #     idepix["conditions/geometry/x"] = x_da
    #     idepix["conditions/geometry/crs"] = crs_da
    #
    #     return geographic_coords_dt

    def _compute_idepix(
            self,
            l1c: xr.DataTree,
            idepix: xr.DataTree,
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int] = None,
    ):
        """Adds quality flags bands in target resolution"""
        band_resolution = self.resolutions[self.bands[0]]
        band_chunksize = chunksize_in_meters // band_resolution
        flag_band_name = "quality/mask/pixel_classif_flag"

        idepix_flags = IdepixClassificationPUPU().run(
            {'l1c': l1c},
        )["pixel_classif_flag"]

        idepix["quality/mask/pixel_classif_flag"] = self._with_encoding(xr.DataArray(
            idepix_flags["pixel_classif_flag"]
        ))


    def _with_encoding(self, var: xr.DataArray) -> xr.DataArray:
        var.encoding = {
            "zlib": True,
            "complevel": 5,
            "chunksizes": var.data.chunksize,
        }
        return var
