from typing import Any, Optional, Literal, Tuple, Dict, Iterable, List

import dask.array as da
import numpy as np
import xarray as xr
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType
from pyproj import Transformer, CRS

from l2w_v1.exceptions.errors import MyError
from sen2water.msiresampling.ancillaryinterpolation import AncillaryInterpolation
from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.constants import MsiConstantsReengineering
from sen2water.msiresampling.downsampling_pu import DownsamplingPU
from sen2water.msiresampling.geocoordinates import GeoCoordinates
from sen2water.msiresampling.meanangles import MeanAngles
from sen2water.msiresampling.tpinterpolation import TpInterpolation
from sen2water.msiresampling.upsampling import Upsampling
from sen2water.msiresampling.upsampling_pu import UpsamplingPU


class ResamplingMainPU(EOProcessingUnit):
    PROCESSOR_NAME = "MsiResampling2"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    INPUT_IDENTIFIER: str = "l1c"

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
            merge_flags: bool = False,
            **kwargs
    ) -> MappingDataType:
        if len(inputs) < 1:
            raise MyError("Missing mandatory input product 'l1c'")

        l1c: xr.DataTree = inputs[self.INPUT_IDENTIFIER]

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {
            "y": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["x"]
        }
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)
        resampled = xr.DataTree(name="resampled")
        # resampled["measurements/reflectance/resampled/x"]

        # copy 1-d metric coordinate bands

        y = xr.DataArray(
            l1c[f"measurements/reflectance/r{resolution}m/y"].data,
            dims=['y'],
            attrs={"long_name": "y coordinate of projection",
                   "standard_name": "projection_y_coordinate",
                   "units": "m"}
        )
        x = xr.DataArray(
            l1c[f"measurements/reflectance/r{resolution}m/x"].data,
            dims=['x'],
            attrs={"long_name": "x coordinate of projection",
                   "standard_name": "projection_x_coordinate",
                   "units": "m"}
        )

        resampled["measurements/reflectance/resampled/y"] = y
        resampled["measurements/reflectance/resampled/x"] = x

        yy, xx = self.construct_coordinate_arrays(
            da.from_array(l1c[f"measurements/reflectance/r{resolution}m/y"].data, chunks=chunks[1]),
            da.from_array(l1c[f"measurements/reflectance/r{resolution}m/x"].data, chunks=chunks[0]),
            chunks=chunks,
        )
        crs = CRS.from_authority(*l1c.attrs["stac_discovery"]["properties"]["proj:code"].split(":"))
        proj_transform = l1c.attrs["stac_discovery"]["properties"]["proj:transform"]
        # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
        # GT(1) w-e pixel resolution / pixel width.
        # GT(2) row rotation (typically zero).
        # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
        # GT(4) column rotation (typically zero).
        # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
        resampled["measurements/reflectance/resampled/crs"] = xr.DataArray(
            0, attrs={
                "crs_wkt": crs.to_wkt(),
                "GeoTransform": " ".join(str(proj_transform[i]) for i in [2, 0, 1, 5, 3, 4])
            }
        )
        transformer = Transformer.from_crs(crs, CRS("EPSG:4326"))
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
        del xx, yy, transformer, lat_lon

        resampled["measurements/reflectance/resampled/latitude"] = xr.DataArray(
            lat_data,
            dims=['y', 'x'],
            attrs={"standard_name": "latitude",
                   "units": "degrees_north"}
        )
        resampled["measurements/reflectance/resampled/longitude"] = xr.DataArray(
            lon_data,
            dims=['y', 'x'],
            attrs={"standard_name": "longitude",
                   "units": "degrees_east"}
        )

        resampled["conditions/mask/detector_footprint/resampled/y"] = resampled["/measurements/reflectance/resampled/y"]
        resampled["conditions/mask/detector_footprint/resampled/x"] = resampled["/measurements/reflectance/resampled/x"]
        resampled["conditions/mask/detector_footprint/resampled/crs"] = resampled[
            "/measurements/reflectance/resampled/crs"]
        if ancillary:
            resampled["conditions/meteorology/resampled/y"] = resampled["/measurements/reflectance/resampled/y"]
            resampled["conditions/meteorology/resampled/x"] = resampled["/measurements/reflectance/resampled/x"]
            resampled["conditions/meteorology/resampled/crs"] = resampled["/measurements/reflectance/resampled/crs"]
        resampled["quality/mask/resampled/y"] = resampled["/measurements/reflectance/resampled/y"]
        resampled["quality/mask/resampled/x"] = resampled["/measurements/reflectance/resampled/x"]
        resampled["quality/mask/resampled/crs"] = resampled["/measurements/reflectance/resampled/crs"]
        resampled["conditions/mask/l1c_classification/resampled/y"] = resampled["/measurements/reflectance/resampled/y"]
        resampled["conditions/mask/l1c_classification/resampled/x"] = resampled["/measurements/reflectance/resampled/x"]
        resampled["conditions/mask/l1c_classification/resampled/crs"] = resampled[
            "/measurements/reflectance/resampled/crs"]
        resampled["conditions/geometry/resampled/y"] = resampled["/measurements/reflectance/resampled/y"]
        resampled["conditions/geometry/resampled/x"] = resampled["/measurements/reflectance/resampled/x"]
        resampled["conditions/geometry/resampled/crs"] = resampled["/measurements/reflectance/resampled/crs"]

        resampled["conditions/meteorology/ecmwf"] = l1c["conditions/meteorology/ecmwf"].copy(inherit=False, deep=False)
        resampled["conditions/meteorology/cams"] = l1c["conditions/meteorology/cams"].copy(inherit=False, deep=False)

        # Removes mismatched encoding var when reading from zarr v2 and writing to zarr v3
        for node in resampled["conditions/meteorology"].subtree:
            for var in node.variables.values():
                # Alternatively: set encoding to a good value
                var.encoding.clear()

        # select detector per target pixel (with_detfoo_filter) or per target pixel per band

        self._resample_detectors(resolution, dims, l1c, resampled, chunksize_in_meters=chunksize_in_meters,
                                 with_detfoo_filter=with_detfoo_filter)

        # resample reflectance bands B1 .. B12

        input_data_with_target_resolution = self._resample_reflectance(
            downsampling,
            upsampling,
            resolution,
            overlap_depth,
            dims,
            l1c,
            resampled,
            chunksize_in_meters=chunksize_in_meters,
            with_detfoo_filter=with_detfoo_filter
        )

        # resample flag bands B_xxx_B1 .. B_zzz_B12

        self._resample_flags(flagdownsampling, resolution, chunksize_in_meters=chunksize_in_meters, dims=dims, l1c=l1c,
                             resampled=resampled)

        self._resample_cloud_ice(flagdownsampling, resolution, chunksize_in_meters=chunksize_in_meters, dims=dims,
                                 l1c=l1c, resampled=resampled)

        if ancillary:
            self._resample_ancillary(
                ancillary,
                l1c,
                lat_data,
                lon_data, dims=dims, resampled=resampled,
            )
        del lat_data, lon_data

        self._resample_sun_angles(
            resolution,
            chunksize_in_meters=chunksize_in_meters,
            dims=dims,
            input_data_with_target_resolution=input_data_with_target_resolution,
            l1c=l1c,
            resampled=resampled
        )

        self._resample_viewing_angles(
            resolution,
            chunksize_in_meters=chunksize_in_meters,
            dims=dims,
            l1c=l1c,
            resampled=resampled,
            with_detfoo_filter=with_detfoo_filter
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

    def _resample_detectors(
            self,
            resolution: int,
            dims: Dict[str, int],
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            chunksize_in_meters: int,
            with_detfoo_filter: bool = True
    ):
        """Adds resampled detector index"""
        if with_detfoo_filter:
            # The result shall be a single detfoo with fill value in the overlapping area without a common detfoo value.
            # Any of the source bands with the target resolution determines the master detector footprint value.
            # If any of the other bands does not have a contribution with this master detfoo value then the master detfoo pixel is set to invalid.
            master_band = None
            master_detfoo = None
            for band in self.bands:
                if self.resolutions[band] == resolution:
                    master_band = band
                    detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{resolution}m/{band}"
                    master_detfoo_data = l1c[detector_footprint_band_name].data.rechunk(
                        chunksize_in_meters // resolution)
                    master_detfoo = xr.DataTree()
                    master_detfoo['master_detfoo'] = xr.DataArray(master_detfoo_data,
                                                                  dims=l1c[detector_footprint_band_name].dims)

                    break
            for band in self.bands:
                if band == master_band:
                    continue
                band_resolution = self.resolutions[band]
                detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
                band_chunksize = chunksize_in_meters // band_resolution
                master_detfoo = DownsamplingPU().run(
                    {'data': l1c, 'target_detfoo': master_detfoo},
                    var_name=detector_footprint_band_name,
                    band_chunksize=band_chunksize,
                    mode="masterdetfoo",
                    factor=factor,
                    dtype=l1c[detector_footprint_band_name].data.dtype,
                    chunks=(chunksize_in_meters // resolution,
                            chunksize_in_meters // resolution)
                )
                if resolution > band_resolution:
                    factor = resolution // band_resolution
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
                        chunks=(band_chunksize * factor, band_chunksize * factor)
                    )
                    master_detfoo['master_detfoo'] = xr.DataArray(da.where(
                        resampled_detector[detector_footprint_band_name].data == master_detfoo['master_detfoo'].data,
                        master_detfoo['master_detfoo'].data, 0), dims=dims)
                else:
                    resampled_detector = l1c
                    master_detfoo['master_detfoo'] = xr.DataArray(da.where(
                        resampled_detector[detector_footprint_band_name].data == master_detfoo['master_detfoo'].data,
                        master_detfoo['master_detfoo'].data, 0), dims=dims)
            resampled["/conditions/mask/detector_footprint/resampled/master_detfoo"] = xr.DataArray(
                master_detfoo,
                dims=dims,
                attrs={"long_name": "master detector footprint",
                       "_FillValue": np.uint8(0),
                       "coordinates": "crs y x"}
            )
        else:
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
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )

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
                        chunks=(band_chunksize * factor, band_chunksize * factor)
                    )
                else:
                    resampled_detector = l1c
                target_band_name = f"/conditions/mask/detector_footprint/resampled/{band}"
                resampled[target_band_name] = xr.DataArray(
                    resampled_detector[detector_footprint_band_name].data,
                    dims=dims,
                    attrs={"long_name": f"detector footprint of {band}",
                           "_FillValue": np.uint8(0),
                           "coordinates": "crs y x"}
                )

    def _resample_reflectance(
            self,
            downsampling: Literal[
                "detectormean", "mean", "median", "min", "max", "first",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
            ],
            upsampling: Literal["nearest", "bilinear", "bicubic"],
            resolution: int,
            overlap_depth: int,
            dims: Dict[str, int],
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            chunksize_in_meters: int,
            with_detfoo_filter: bool,
    ) -> da.Array:
        """Adds resampled reflectance bands, returns one input band with target resolution as dummy"""
        input_data_with_target_resolution: Optional[da.Array] = None

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
                    target_tree['target_detector_index'] = xr.DataArray(
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
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )
                else:
                    resampled_band = DownsamplingPU().run(
                        {'data': l1c},
                        var_name=band_name,
                        band_chunksize=band_chunksize,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=l1c[band_name].data.dtype,
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )
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
                    chunks=(band_chunksize * factor, band_chunksize * factor)
                )
            else:
                resampled_band = l1c
                if input_data_with_target_resolution is None:
                    input_data_with_target_resolution = xr.DataTree()
                    input_data_with_target_resolution['data'] = xr.DataArray(l1c[band_name].data, dims=dims)

            attrs = resampled_band[band_name].attrs
            attrs["coordinates"] = "crs y x latitude longitude"
            resampled[f"measurements/reflectance/resampled/{band}"] = xr.DataArray(
                resampled_band[band_name].data,
                dims=dims,
                attrs=attrs
            )
        return input_data_with_target_resolution

    def _resample_ancillary(
            self,
            ancillary: Iterable[str],
            l1c: xr.DataTree,
            lat_data: da.Array,
            lon_data: da.Array,
            dims: Dict[str, int],
            resampled: xr.DataTree
    ):
        for anc_band_name in ancillary:
            anc_lat = l1c["conditions/meteorology/ecmwf/latitude"]
            anc_lon = l1c["conditions/meteorology/ecmwf/longitude"]
            anc_band = l1c[anc_band_name]
            anc_data = AncillaryInterpolation().apply(
                lat_data,
                lon_data,
                anc_lat=anc_lat.values,
                anc_lon=anc_lon.values,
                anc_data=anc_band.values,
                variable=anc_band,
                dtype=np.float32,
            )
            attrs = anc_band.attrs.copy()
            attrs["coordinates"] = "crs y x"
            resampled[f"conditions/meteorology/resampled/{anc_band.name}"] = xr.DataArray(anc_data, dims=dims,
                                                                                          attrs=attrs)

    def _resample_flags(
            self,
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            l1c: xr.DataTree,
            resampled: xr.DataTree,
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
                    chunks=(chunksize_in_meters // resolution,
                            chunksize_in_meters // resolution)
                )
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
                    chunks=(band_chunksize * factor, band_chunksize * factor)
                )
            else:
                resampled_flags = l1c
            attrs = l1c[flag_band_name].attrs.copy()
            attrs["coordinates"] = "crs y x"
            resampled[f"quality/mask/resampled/{band}"] = xr.DataArray(
                resampled_flags[flag_band_name].data,
                dims=dims,
                attrs=attrs
            )

    def _resample_cloud_ice(
            self,
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            l1c: xr.DataTree,
            resampled: xr.DataTree
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
                chunks=(band_chunksize * factor, band_chunksize * factor)
            )
        else:
            resampled_flags = l1c
        attrs = l1c[band_name].attrs.copy()
        attrs["coordinates"] = "crs y x"

        resampled[f"conditions/mask/l1c_classification/resampled/b00"] = xr.DataArray(
            resampled_flags[band_name].data,
            dims=dims,
            attrs=attrs
        )

    def _resample_sun_angles(
            self,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            input_data_with_target_resolution: da.Array,
            l1c: xr.DataTree,
            resampled: xr.DataTree
    ):
        """Adds sun angles resampled from TP grid in target resolution"""

        angle_indices = {
            0: ("sza", "Sun Zenith Angle"),
            1: ("saa", "Sun Azimuth Angle"),
        }
        band = "conditions/geometry/sun_angles"
        for angle_idx, (angle_short_name, angle_long_name) in angle_indices.items():
            angle_data = l1c[band].data[angle_idx]
            band_chunksize = chunksize_in_meters // resolution
            resampled_angles = TpInterpolation().apply(
                input_data_with_target_resolution,
                tp_data=angle_data,
                resolution=resolution,
                tp_resolution=5000,
                image_shape=input_data_with_target_resolution.shape,
                image_chunksize=(band_chunksize, band_chunksize),
                dtype=angle_data.dtype,
            )
            resampled[f"conditions/geometry/resampled/{angle_short_name}"] = xr.DataArray(
                resampled_angles,
                dims=dims,
                attrs={"long_name": angle_long_name,
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x"}
            )

    def _resample_viewing_angles(
            self,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            with_detfoo_filter: bool = True,
    ):
        """Resamples viewing angles per detector and adds viewing angles per band"""
        vza_accu = []
        vaa_accu = []

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
        target_chunksize = chunksize_in_meters // resolution
        detectors = np.array(
            [int(d.removeprefix("d")) for d in view_angle_band.coords["detector"].values]
        )

        for band in self.bands:
            detector_footprint_band_name = "conditions/mask/detector_footprint/resampled/" + (
                "master_detfoo" if with_detfoo_filter else band)
            detector_footprint = resampled[detector_footprint_band_name].data
            extended_vza, extended_vaa = AnglesInterpolation().expand_angles_per_detector(
                view_angle_band.isel(angle=angle_indices["vza"]).sel(band=band).values,
                view_angle_band.isel(angle=angle_indices["vaa"]).sel(band=band).values,
            )
            resampled_vza = AnglesInterpolation().apply(
                detector_footprint,
                detectors=detectors,
                detector_angles=extended_vza,
                resolution=resolution,
                angles_resolution=5000,
                band=band,
                image_chunksize=(target_chunksize, target_chunksize),
                dtype=np.float32,  # TBC
            )
            resampled_vaa = AnglesInterpolation().apply(
                detector_footprint,
                detectors=detectors,
                detector_angles=extended_vaa,
                resolution=resolution,
                angles_resolution=5000,
                band=band,
                is_azimuth_angle=True,
                image_chunksize=(target_chunksize, target_chunksize),
                dtype=np.float32,  # TBC
            )
            # Shape of the view_angle_band
            # "coordinates": [
            #   "detector",
            #   "angle",
            #   "y",
            #   "x",
            #   "band"
            # ],
            resampled[f"conditions/geometry/resampled/vza_{band}"] = xr.DataArray(
                resampled_vza,
                dims=dims,
                attrs={"long_name": f"Viewing incidence zenith angle for band {band}",
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x"}
            )
            resampled[f"conditions/geometry/resampled/vaa_{band}"] = xr.DataArray(
                resampled_vaa,
                dims=dims,
                attrs={"long_name": f"Viewing incidence azimuth angle for band {band}",
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x"}
            )
            vza_accu.append(resampled_vza)
            vaa_accu.append(resampled_vaa)

        vza_mean = MeanAngles().apply(*vza_accu, dtype=np.float32)
        vaa_mean = MeanAngles().apply(*vaa_accu, is_azimuth_angle=True, dtype=np.float32)

        resampled[f"conditions/geometry/resampled/vza_mean"] = xr.DataArray(
            vza_mean,
            dims=dims,
            attrs={"long_name": "Mean viewing incidence zenith angle",
                   "units": "degrees",
                   "_FillValue": np.nan,
                   "coordinates": "crs y x"}
        )
        resampled[f"conditions/geometry/resampled/vaa_mean"] = xr.DataArray(
            vaa_mean,
            dims=dims,
            attrs={"long_name": "Mean viewing incidence azimuth angle",
                   "units": "degrees",
                   "_FillValue": np.nan,
                   "coordinates": "crs y x"}
        )

    # @dask.delayed
    def construct_coordinate_arrays(self, y, x, chunks: Tuple[int, int]):
        xx = da.broadcast_to(x, (y.shape[0], x.shape[0]), chunks=chunks)
        yy = da.transpose(da.broadcast_to(y, (y.shape[0], x.shape[0]), chunks=chunks[::-1]))
        return yy, xx
