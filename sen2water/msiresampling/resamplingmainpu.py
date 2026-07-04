from typing import Any, Optional, Literal, Tuple, Dict, Iterable, List
import numpy as np
import xarray as xr
from dask import array as da
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType
from pyproj import Transformer, CRS
from sen2water.msiresampling.ancillaryinterpolation_pu import AncillaryInterpolationPU
from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.anglesinterpolation_pu import AnglesInterpolationPU
from sen2water.msiresampling.constants import MsiConstantsReengineering
from sen2water.msiresampling.downsampling_pu import DownsamplingPU
from sen2water.msiresampling.geocoordinates_pu import GeoCoordinatesPU
from sen2water.msiresampling.meanangles_pu import MeanAnglesPU
from sen2water.msiresampling.tpinterpolation_pu import TpInterpolationPU
from sen2water.msiresampling.upsampling_pu import UpsamplingPU


class ResamplingMainPU(EOProcessingUnit):
    PROCESSOR_NAME = "MsiResampling2"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    bands = MsiConstantsReengineering.bands
    resolutions = MsiConstantsReengineering.resolutions

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
            "y": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["x"]
        }
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)

        # create result accumulator

        resampled = xr.DataTree(name="resampled")

        # add geo-coding with 1-d metric coordinates, CRS, 2-d geographic coordinates

        geographic_coords_dt = self._add_geocoding(
            l1c,
            resampled,
            ancillary,
            resolution,
            dims,
            chunks,
        )

        # select detector per target pixel (with_detfoo_filter) or per target pixel per band

        if with_detfoo_filter:
            self._resample_detectors_with_detfoo_filter(
                l1c, resampled,
                resolution, chunksize_in_meters, dims, chunks,
            )
        else:
            self._resample_detectors(
                l1c, resampled,
                resolution, chunksize_in_meters, dims, chunks,
            )

        # resample reflectance bands b01 .. b12

        input_data_with_target_resolution = self._resample_reflectance(
            l1c, resampled,
            downsampling, upsampling, with_detfoo_filter, overlap_depth,
            resolution, chunksize_in_meters, dims, chunks,
        )

        # resample flag bands B_xxx_B1 .. B_zzz_B12

        self._resample_flags(
            l1c, resampled,
            flagdownsampling,
            resolution, chunksize_in_meters,  dims, chunks
        )

        self._resample_cloud_ice(
            l1c, resampled,
            flagdownsampling,
            resolution, chunksize_in_meters, dims, chunks
        )

        # copy coarse grained meteo data

        resampled["conditions/meteorology/ecmwf"] = l1c["conditions/meteorology/ecmwf"].copy(inherit=False, deep=False)
        resampled["conditions/meteorology/cams"] = l1c["conditions/meteorology/cams"].copy(inherit=False, deep=False)

        # Removes mismatched encoding var when reading from zarr v2 and writing to zarr v3
        for node in resampled["conditions/meteorology"].subtree:
            for var in node.variables.values():
                var.encoding.clear()

        # resample meteo data

        if ancillary:
            self._resample_ancillary(
                l1c, resampled,
                ancillary, geographic_coords_dt,
                dims, chunks,
            )

        # resample angles

        self._resample_sun_angles(
            l1c, input_data_with_target_resolution, resampled,
            resolution, chunksize_in_meters, dims, chunks,
        )

        self._resample_viewing_angles(
            l1c, resampled,
            with_detfoo_filter,
            resolution, chunksize_in_meters, dims, chunks,
        )

        properties = l1c.attrs["stac_discovery"]["properties"]
        resampled.attrs = {
            "other_metadata": {
                "product_type": "S2_MSI_Level-1C",
                "platform": properties["platform"],
                "start_date": properties["start_datetime"],
                "stop_date": properties["end_datetime"],
            }
        }

        return {
            "resampled": resampled
        }

    def _add_geocoding(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            ancillary: List[str]|None,
            resolution: int,
            dims: dict[str, int],
            chunks: tuple[int, int],
    ) -> xr.DataTree:

        epsg_code = l1c.attrs["stac_discovery"]["properties"]["proj:code"]
        # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
        # GT(1) w-e pixel resolution / pixel width.
        # GT(2) row rotation (typically zero).
        # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
        # GT(4) column rotation (typically zero).
        # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
        proj_transform = l1c.attrs["stac_discovery"]["properties"]["proj:transform"]

        geo_transform = " ".join(str(proj_transform[i]) for i in [2, 0, 1, 5, 3, 4])
        crs = CRS.from_authority(*epsg_code.split(":"))
        transformer = Transformer.from_crs(crs, CRS("EPSG:4326"))

        # The metric coordinates are numpy arrays, they are written as a single chunk.
        y_1d = l1c[f"measurements/reflectance/r{resolution}m/y"].data
        x_1d = l1c[f"measurements/reflectance/r{resolution}m/x"].data
        # Here we construct dask arrays to make the 2-d lat and lon chunked.
        y_data = da.from_array(y_1d, chunks=chunks[1])
        x_data = da.from_array(x_1d, chunks=chunks[0])
        y_grid = da.transpose(da.broadcast_to(y_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks[::-1]))
        x_grid = da.broadcast_to(x_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks)

        metric_coords_dt = xr.DataTree()
        metric_coords_dt["y"] = xr.DataArray(y_grid, dims=dims, )
        metric_coords_dt["x"] = xr.DataArray(x_grid, dims=dims, )

        geographic_coords_dt = GeoCoordinatesPU().run(
            {"metric_coords": metric_coords_dt},
            transformer=transformer,
            dtype=np.float64,
            dims=dims,
            chunks=chunks,
        )

        y_da = xr.DataArray(
            y_1d,
            dims=["y"],
            attrs={
                "long_name": "y coordinate of projection",
                "standard_name": "projection_y_coordinate",
                "units": "m",
            },
        )
        x_da = xr.DataArray(
            x_1d,
            dims=["x"],
            attrs={
                "long_name": "x coordinate of projection",
                "standard_name": "projection_x_coordinate",
                "units": "m",
            },
        )
        crs_da = xr.DataArray(
            0,
            attrs={
                "crs_wkt": crs.to_wkt(),
                "GeoTransform": geo_transform
            }
        )
        lat_da = geographic_coords_dt["geographic_coords"]["latitude"]
        lon_da = geographic_coords_dt["geographic_coords"]["longitude"]

        resampled["measurements/reflectance/resampled/y"] = y_da
        resampled["measurements/reflectance/resampled/x"] = x_da
        resampled["measurements/reflectance/resampled/crs"] = crs_da
        resampled["measurements/reflectance/resampled/latitude"] = self._with_encoding(lat_da)
        resampled["measurements/reflectance/resampled/longitude"] = self._with_encoding(lon_da)
        resampled["conditions/mask/detector_footprint/resampled/y"] = y_da
        resampled["conditions/mask/detector_footprint/resampled/x"] = x_da
        resampled["conditions/mask/detector_footprint/resampled/crs"] = crs_da
        if ancillary:
            resampled["conditions/meteorology/resampled/y"] = y_da
            resampled["conditions/meteorology/resampled/x"] = x_da
            resampled["conditions/meteorology/resampled/crs"] = crs_da
        resampled["conditions/mask/l1c_classification/resampled/y"] = y_da
        resampled["conditions/mask/l1c_classification/resampled/x"] = x_da
        resampled["conditions/mask/l1c_classification/resampled/crs"] = crs_da
        resampled["conditions/geometry/resampled/y"] = y_da
        resampled["conditions/geometry/resampled/x"] = x_da
        resampled["conditions/geometry/resampled/crs"] = crs_da
        resampled["quality/mask/resampled/y"] = y_da
        resampled["quality/mask/resampled/x"] = x_da
        resampled["quality/mask/resampled/crs"] = crs_da

        return geographic_coords_dt

    def _resample_detectors_with_detfoo_filter(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: tuple[int, int],
    ):
        """Adds resampled detector index"""
        # The result shall be a single detfoo with fill value in the overlapping area without a common detfoo value.
        # Any of the source bands with the target resolution determines the master detector footprint value.
        # If any of the other bands does not have a contribution with this master detfoo value then the master detfoo pixel is set to invalid.
        master_band, master_detfoo = self._initialise_master_detfoo(l1c, resolution, chunks)

        for band in self.bands:
            if band == master_band:
                continue
            band_resolution = self.resolutions[band]
            detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
            band_chunksize = chunksize_in_meters // band_resolution

            if resolution > band_resolution:
                factor = resolution // band_resolution
                master_detfoo = DownsamplingPU().run(
                    {'data': l1c, 'target_detfoo': master_detfoo},
                    var_name=detector_footprint_band_name,
                    band_chunksize=band_chunksize,
                    mode="masterdetfoo",
                    factor=factor,
                    dtype=l1c[detector_footprint_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]

            elif resolution < band_resolution:
                factor = band_resolution // resolution
                resampled_detector = UpsamplingPU().run(
                    {'data': l1c},
                    var_name=detector_footprint_band_name,
                    band_chunksize=band_chunksize,
                    src_image_shape=l1c[detector_footprint_band_name].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    mode="nearest",
                    factor=(factor, factor),
                    depth=0,
                    dtype=l1c[detector_footprint_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]
                # pixels may be not in line with the master
                master_detfoo['target_detfoo'] = xr.DataArray(
                    da.where(
                        resampled_detector[detector_footprint_band_name].data == master_detfoo['target_detfoo'].data,
                        master_detfoo['target_detfoo'].data,
                        0),
                    dims=dims,
                    attrs=resampled_detector[detector_footprint_band_name].attrs | {"coordinates": "y x"},  # TODO check attrs of input
                )

            else:
                resampled_detector_data = l1c[detector_footprint_band_name].data.rechunk(chunks)
                # pixels may be not in line with the master
                master_detfoo['target_detfoo'] = xr.DataArray(
                    da.where(
                        resampled_detector_data == master_detfoo['target_detfoo'].data,
                        master_detfoo['target_detfoo'].data,
                        0),
                    dims=dims,
                    attrs=l1c[detector_footprint_band_name].attrs | {
                        "coordinates": "y x",
                        "grid_mapping": "crs",
                    },  # TODO check attrs of input
                )

        resampled["/conditions/mask/detector_footprint/resampled/master_detfoo"] = (
            self._with_encoding(master_detfoo["target_detfoo"])
        )

    def _initialise_master_detfoo(
            self, l1c: xr.DataTree, resolution: int, chunks: tuple[int, int]
    ) -> tuple[str, xr.DataTree]:
        master_band = None
        master_detfoo = None
        for band in self.bands:
            if self.resolutions[band] == resolution:
                master_band = band
                detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{resolution}m/{band}"
                master_detfoo_data = l1c[detector_footprint_band_name].data.rechunk(chunks)
                master_detfoo = xr.DataTree()
                master_detfoo['target_detfoo'] = xr.DataArray(
                    master_detfoo_data,
                    dims=l1c[detector_footprint_band_name].dims,
                    attrs={"long_name": "master detector footprint",
                           "_FillValue": np.uint8(0),
                           "coordinates": "y x",
                           "grid_mapping": "crs",
                    }
                )
                break
        if master_detfoo is None:
            raise KeyError(f"Missing input band with resolution {resolution} for master detfoo.")
        return master_band, master_detfoo

    def _resample_detectors(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: tuple[int, int],
    ):
        """Adds resampled detector index"""
        for band in self.bands:
            band_resolution = self.resolutions[band]
            detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
            band_chunksize = chunksize_in_meters // band_resolution

            if resolution > band_resolution:
                factor = resolution // band_resolution
                resampled_detector = DownsamplingPU().run(
                    {'data': l1c},
                    var_name=detector_footprint_band_name,
                    band_chunksize=band_chunksize,
                    mode="majority",
                    factor=factor,
                    dtype=l1c[detector_footprint_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]

            elif resolution < band_resolution:
                factor = band_resolution // resolution
                resampled_detector = UpsamplingPU().run(
                    {'data': l1c},
                    var_name=detector_footprint_band_name,
                    band_chunksize=band_chunksize,
                    mode="nearest",
                    factor=(factor, factor),
                    src_image_shape=l1c[detector_footprint_band_name].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=l1c[detector_footprint_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]
            else:
                resampled_detector = xr.DataTree()
                resampled_detector[detector_footprint_band_name] = xr.DataArray(
                    l1c[detector_footprint_band_name].data.rechunk(chunks),
                    dims=dims,
                    attrs=l1c[detector_footprint_band_name].attrs,
                )

            resampled_detector_da = xr.DataArray(
                resampled_detector[detector_footprint_band_name].data,
                dims=dims,
                attrs={
                    "long_name": f"detector footprint of {band}",
                    "_FillValue": np.uint8(0),
                    "coordinates": "y x",
                    "grid_mapping": "crs",}
            )
            resampled_detector_da.encoding = {
                    "zlib": True,
                    "complevel": 5,
                    "chunksizes": resampled_detector_da.data.chunksize,
            }
            resampled[f"/conditions/mask/detector_footprint/resampled/{band}"] = resampled_detector_da

    def _resample_reflectance(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            downsampling: Literal[
                "detectormean",
                "mean",
                "median",
                "min",
                "max",
                "first",
                "flagand",
                "flagor",
                "flagmedianand",
                "flagmedianor",
            ],
            upsampling: Literal["nearest", "bilinear", "bicubic"],
            with_detfoo_filter: bool,
            overlap_depth: int,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ) -> xr.DataTree:
        """Adds resampled reflectance bands, returns one input band with target resolution as dummy"""
        input_data_with_target_resolution: Optional[xr.DataTree] = None

        for band in self.bands:
            band_resolution = self.resolutions[band]
            band_chunksize = chunksize_in_meters // band_resolution
            band_name = f"measurements/reflectance/r{band_resolution}m/{band}"
            factor = resolution // band_resolution
            if resolution > band_resolution:

                if downsampling == 'detectormean':
                    detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
                    detector_footprint_band_data = l1c[detector_footprint_band_name].data.rechunk(
                        chunksize_in_meters // band_resolution)

                    target_tree = xr.DataTree()
                    target_tree['target_detfoo'] = xr.DataArray(
                        resampled["conditions/mask/detector_footprint/resampled/" + (
                            "master_detfoo" if with_detfoo_filter else band)].data,
                        dims=l1c[detector_footprint_band_name].dims
                    )

                    detector_tree = xr.DataTree()
                    detector_tree['detfoo'] = xr.DataArray(detector_footprint_band_data,
                                                           dims=l1c[detector_footprint_band_name].dims)

                    resampled_band = DownsamplingPU().run(
                        {'data': l1c,
                         'target_detfoo': target_tree,
                         'detfoo': detector_tree},
                        var_name=band_name,
                        band_chunksize=band_chunksize,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=l1c[band_name].data.dtype,
                        chunks=chunks
                    )["resampled"]
                else:
                    resampled_band = DownsamplingPU().run(
                        {'data': l1c},
                        var_name=band_name,
                        band_chunksize=band_chunksize,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=l1c[band_name].data.dtype,
                        chunks=chunks
                    )["resampled"]
            elif resolution < self.resolutions[band]:
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_band = UpsamplingPU().run(
                    {'data': l1c},
                    var_name=band_name,
                    band_chunksize=band_chunksize,
                    mode=upsampling,
                    factor=(factor, factor),
                    src_image_shape=l1c[band_name].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    is_reflectance=True,
                    depth=overlap_depth,
                    trim=False,
                    dtype=l1c[band_name].data.dtype,
                    chunks=chunks
                )["resampled"]
            else:
                resampled_band = xr.DataTree()
                resampled_band[band_name] = xr.DataArray(
                    l1c[band_name].data.rechunk(chunks),
                    dims=dims,
                    attrs=l1c[band_name].attrs,
                )
                if input_data_with_target_resolution is None:
                    input_data_with_target_resolution = xr.DataTree()
                    input_data_with_target_resolution["data"] = resampled_band[band_name]

            resampled_band_da = xr.DataArray(
                resampled_band[band_name].data,
                dims=dims,
                attrs=resampled_band[band_name].attrs | {
                    "coordinates": "y x latitude longitude",
                    "grid_mapping": "crs",
                },
            )
            resampled[f"measurements/reflectance/resampled/{band}"] = self._with_encoding(resampled_band_da)
        return input_data_with_target_resolution

    def _resample_ancillary(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            ancillary: Iterable[str],
            geographic_coords: xr.DataTree,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ):
        anc_lat = l1c["conditions/meteorology/ecmwf/latitude"].values
        anc_lon = l1c["conditions/meteorology/ecmwf/longitude"].values
        for anc_band_name in ancillary:
            anc_result = AncillaryInterpolationPU().run(
                {"l1c": l1c, "geographic_coords": geographic_coords["geographic_coords"]},
                anc_lat=anc_lat,
                anc_lon=anc_lon,
                anc_band_name=anc_band_name,
                dtype=np.float32,
                chunks=chunks,
            )
            resampled[f"conditions/meteorology/resampled/{anc_band_name}"] = self._with_encoding(
                anc_result["ancillary"][anc_band_name]
            )

    def _resample_flags(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
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
        for band in self.bands:
            band_resolution = self.resolutions[band]
            band_chunksize = chunksize_in_meters // band_resolution
            flag_band_name = f"quality/mask/r{band_resolution}m/{band}"

            if resolution > band_resolution:
                factor = resolution // band_resolution
                resampled_flags = DownsamplingPU().run(
                    {'data': l1c},
                    var_name=flag_band_name,
                    band_chunksize=band_chunksize,
                    mode=flagdownsampling,
                    factor=factor,
                    dtype=l1c[flag_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]
            elif resolution < band_resolution:
                factor = band_resolution // resolution
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_flags = UpsamplingPU().run(
                    {'data': l1c},
                    var_name=flag_band_name,
                    band_chunksize=band_chunksize,
                    mode='nearest',
                    factor=(factor, factor),
                    src_image_shape=l1c[flag_band_name].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=l1c[flag_band_name].data.dtype,
                    chunks=chunks,
                )["resampled"]
            else:
                resampled_flags = xr.DataTree()
                resampled_flags[flag_band_name] = xr.DataArray(
                    l1c[flag_band_name].data.rechunk(chunks),
                    dims=dims,
                    attrs=l1c[flag_band_name].attrs
                )
            resampled[f"quality/mask/resampled/{band}"] = self._with_encoding(xr.DataArray(
                resampled_flags[flag_band_name].data,
                dims=dims,
                attrs=l1c[flag_band_name].attrs.copy() | {
                    "coordinates": "y x",
                    "gird_mapping": "crs",
                }
            ))

    def _resample_cloud_ice(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
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
        """Adds cloud and ice flags in target resolution"""
        band_name = "conditions/mask/l1c_classification/r60m/b00"
        band_chunksize = chunksize_in_meters // 60
        if resolution < 60:
            factor = 60 // resolution
            resampled_flags = UpsamplingPU().run(
                {'data': l1c},
                var_name=band_name,
                band_chunksize=band_chunksize,
                mode='nearest',
                factor=(factor, factor),
                src_image_shape=l1c[band_name].data.shape,
                src_image_chunksize=(band_chunksize, band_chunksize),
                depth=0,
                dtype=l1c[band_name].data.dtype,
                chunks=chunks
            )["resampled"]
        else:
            resampled_flags = xr.DataTree()
            resampled_flags[band_name] = xr.DataArray(
                l1c[band_name].data.rechunk(chunks),
                dims=dims,
                attrs=l1c[band_name].attrs
            )

        resampled_flags_da = xr.DataArray(
            resampled_flags[band_name].data,
            dims=dims,
            attrs=l1c[band_name].attrs.copy() | {
                "coordinates": "y x",
                "grid_mapping": "crs",
            }
        )
        resampled[f"conditions/mask/l1c_classification/resampled/b00"] = self._with_encoding(resampled_flags_da)


    def _resample_sun_angles(
            self,
            l1c: xr.DataTree,
            input_data_with_target_resolution: xr.DataTree,
            resampled: xr.DataTree,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ):
        """Adds sun angles resampled from TP grid in target resolution"""

        angle_indices = {
            0: ("sza", "Sun Zenith Angle"),
            1: ("saa", "Sun Azimuth Angle"),
        }
        for angle_idx, (angle_short_name, angle_long_name) in angle_indices.items():
            band_name = "conditions/geometry/sun_angles"
            resampled_angles = TpInterpolationPU().run(
                { "l1c": l1c, "dummy_grid": input_data_with_target_resolution },
                band_name=band_name,
                band_idx=angle_idx,
                band_long_name=angle_long_name,
                resolution=resolution,
                tp_resolution=5000,
                image_chunksize=chunks,
                chunks=chunks,
            )
            resampled[f"conditions/geometry/resampled/{angle_short_name}"] = self._with_encoding(
                resampled_angles["interpolated"][band_name]
            )


    def _resample_viewing_angles(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            with_detfoo_filter: bool,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ):
        """Resamples viewing angles per detector and adds viewing angles per band"""
        vza_accu = dict()
        vaa_accu = dict()

        view_angle_band_name = "conditions/geometry/viewing_incidence_angles"
        view_angle_band = l1c[view_angle_band_name]

        # Shape of the view_angle_band
        # "coordinates": [
        #   "detector",
        #   "angle",
        #   "y",
        #   "x",
        #   "band"
        # ],
        angle_indices = {
            "vza": 0,
            "vaa": 1
        }
        detectors = np.array(
            [int(d.removeprefix("d")) for d in view_angle_band.coords["detector"].values]
        )

        for band in self.bands:
            detector_footprint_band_name = "conditions/mask/detector_footprint/resampled/" + (
                "master_detfoo" if with_detfoo_filter else band)
            extended_vza, extended_vaa = AnglesInterpolation().expand_angles_per_detector(
                view_angle_band.isel(angle=angle_indices["vza"]).sel(band=band).values,
                view_angle_band.isel(angle=angle_indices["vaa"]).sel(band=band).values,
            )
            resampled_vza = AnglesInterpolationPU().run(
                { "detector_footprint": resampled },
                detector_footprint_band_name=detector_footprint_band_name,
                detectors=detectors,
                detector_angles=extended_vza,
                resolution=resolution,
                angles_resolution=5000,
                band_name=band,
                long_name=f"Viewing incidence zenith angle for band {band}",
                image_chunksize=chunks,
                chunks=chunks,
            )
            resampled_vaa = AnglesInterpolationPU().run(
                { "detector_footprint": resampled },
                detector_footprint_band_name=detector_footprint_band_name,
                detectors=detectors,
                detector_angles=extended_vaa,
                resolution=resolution,
                angles_resolution=5000,
                band_name=band,
                long_name=f"Viewing incidence azimuth angle for band {band}",
                is_azimuth_angle=True,
                image_chunksize=chunks,
                chunks=chunks,
            )
            # Shape of the view_angle_band
            # "coordinates": [
            #   "detector",
            #   "angle",
            #   "y",
            #   "x",
            #   "band"
            # ],
            resampled[f"conditions/geometry/resampled/vza_{band}"] = self._with_encoding(
                resampled_vza["interpolated"][band]
            )
            resampled[f"conditions/geometry/resampled/vaa_{band}"] = self._with_encoding(
                resampled_vaa["interpolated"][band]
            )
            vza_accu[band] = resampled_vza["interpolated"]
            vaa_accu[band] = resampled_vaa["interpolated"]

        vza_mean = MeanAnglesPU().run(
            vza_accu,
            band_prefix="conditions/geometry/resampled/vza_",
            long_name="Mean viewing incidence zenith angle",
            dims=dims,
            dtype=np.float32)
        vaa_mean = MeanAnglesPU().run(
            vaa_accu,
            band_prefix="conditions/geometry/resampled/vaa_",
            long_name="Mean viewing incidence azimuth angle",
            dims=dims,
            is_azimuth_angle=True,
            dtype=np.float32)

        resampled["conditions/geometry/resampled/vza_mean"] = self._with_encoding(
                vza_mean["interpolated"]["conditions/geometry/resampled/vza_mean"]
        )
        resampled["conditions/geometry/resampled/vaa_mean"] = self._with_encoding(
                vaa_mean["interpolated"]["conditions/geometry/resampled/vaa_mean"]
        )


    def _with_encoding(self, var: xr.DataArray) -> xr.DataArray:
        var.encoding = {
            "zlib": True,
            "complevel": 5,
            "chunksizes": var.data.chunksize,
        }
        return var
