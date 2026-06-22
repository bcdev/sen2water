from typing import Any, Optional, Literal, Tuple, Dict, Iterable, List

import numpy as np
import xarray as xr
import dask.array as da
from pyproj import Transformer, CRS

from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from l2w_v1.exceptions.errors import MyError
from sen2water.msiresampling.ancillaryinterpolation import AncillaryInterpolation
from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.constants import MsiConstantsReengineering
from sen2water.msiresampling.downsampling import Downsampling
from sen2water.msiresampling.geocoordinates import GeoCoordinates
from sen2water.msiresampling.meanangles import MeanAngles
from sen2water.msiresampling.tpinterpolation import TpInterpolation
from sen2water.msiresampling.upsampling import Upsampling


class ResamplingPU(EOProcessingUnit):

    PROCESSOR_NAME = "MsiResampling"
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
        chunksize_in_meters = 36600,
        downsampling: Literal['detectormean', 'first', 'min', 'max', 'mean', 'median'] = "detectormean",
        flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ] = "flagand",
        upsampling: Literal["nearest", "bilinear", "bicubic"] = "nearest",
        ancillary: List[str],
        with_detfoo_filter: bool = False,
        merge_flags: bool = False,
        **kwargs
    ) -> MappingDataType:
        if len(inputs) < 1:
            raise MyError("Missing mandatory input product 'l1c'")

        l1c: xr.DataTree = inputs[self.INPUT_IDENTIFIER]

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {"y": l1c.sizes[f"y{resolution}"], "x": l1c.sizes[f"x{resolution}"]}
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)
        resampled = xr.DataTree(name="resampled")
        #resampled["measurements/reflectance/resampled/x"]


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
            l1c[f"measurements/reflectance/r{resolution}m/y"],
            l1c[f"measurements/reflectance/r{resolution}m/x"],
            chunks=chunks,
        )
        crs = CRS.from_epsg(l1c.attrs["stac_discovery"]["properties"]["proj:code"])
        proj_transform = l1c.attrs["stac_discovery"]["properties"]["proj:transform"].split(",")
        # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
        # GT(1) w-e pixel resolution / pixel width.
        # GT(2) row rotation (typically zero).
        # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
        # GT(4) column rotation (typically zero).
        # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
        resampled["measurements/reflectance/resampled/crs"] = xr.DataArray(
            0, attrs={
                "crs_wkt": crs.to_wkt(),
                "GeoTransform": " ".join(proj_transform[i] for i in [2, 0, 1, 5, 3, 4])
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

        resampled["measurements/reflectance/resampled/lat"] = xr.DataArray(
            lat_data,
            dims=['y', 'x'],
            attrs={"standard_name": "latitude",
                   "units": "degrees_north"}
        )
        resampled["measurements/reflectance/resampled/lon"] = xr.DataArray(
            lon_data,
            dims=['y', 'x'],
            attrs={"standard_name": "longitude",
                   "units": "degrees_east"}
        )

        resampled["conditions/meteorology/ecmwf"] = l1c["conditions/meteorology/ecmwf"].copy(inherit=False, deep=False)

        # select detector per target pixel (with_detfoo_filter) or per target pixel per band

        self._resample_detectors(resolution, dims, l1c, resampled, chunksize_in_meters=chunksize_in_meters, with_detfoo_filter=with_detfoo_filter)

        # resample reflectance bands B1 .. B12

        input_band_with_target_resolution = self._resample_reflectance(
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

        self._resample_flags(flagdownsampling, resolution, chunksize_in_meters=chunksize_in_meters, dims=dims, l1c=l1c, resampled=resampled)

        self._resample_cloud_ice(flagdownsampling, resolution, chunksize_in_meters=chunksize_in_meters, dims=dims, l1c=l1c, resampled=resampled)

        if ancillary:
            self._resample_ancillary(
                ancillary,
                l1c,
                lat_data,
                lon_data, dims=dims, resampled=resampled,
            )
        del lat_data, lon_data

        # TODO ancillary bands

        self._resample_sun_angles(
            resolution,
            chunksize_in_meters=chunksize_in_meters,
            dims=dims,
            input_band_with_target_resolution=input_band_with_target_resolution,
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

        # TODO check: Copied verbating from ResamplingOperator, only l1c.attrs accesses modified.
        properties = l1c.attrs["stac_discovery"]["properties"]
        resampled.attrs = {
            "conventions": "CF-1.10",
            "TileSize": "610:610",
            "product_type": "S2_MSI_Level-1C",
            "platform": properties["platform"],
            "metadata_profile": "beam",
            "metadata_version": "0.5",
            "auto_grouping": "sun:" \
                             "view:" \
                             "quality:" \
                             "ECMWF:" \
                             "tile:" \
                             "detector_footprint:" \
                             "nodata:" \
                             "partially_corrected_crosstalk:" \
                             "coarse_cloud:" \
                             "snow_and_ice_areas:" \
                             "saturated_l1a:" \
                             "saturated_l1b:" \
                             "defective:ancillary_lost:" \
                             "ancillary_degraded:" \
                             "msi_lost:" \
                             "msi_degraded:" \
                             "saturated_l1a:" \
                             "opaque_clouds:" \
                             "cirrus_clouds:" \
                             "scl:" \
                             "msc:" \
                             "ddv:" \
                             "tile:" \
                             "detector_footprint-B01:detector_footprint-B02:detector_footprint-B03:" \
                             "detector_footprint-B04:detector_footprint-B05:detector_footprint-B06:" \
                             "detector_footprint-B07:detector_footprint-B08:detector_footprint-B8A:" \
                             "detector_footprint-B09:detector_footprint-B10:detector_footprint-B11:" \
                             "detector_footprint-B12:" \
                             "quality_mask",
            "start_date": properties["start_datetime"],
            "stop_date": properties["stop_datetime"],
        }

        return {
            "resampled": resampled
        }


    def _resample_detectors(
        self,
        resolution: int,
        dims: Dict[str,int],
        l1c: xr.DataTree,
        resampled: xr.DataTree,
        chunksize_in_meters: int,
        with_detfoo_filter: bool=True
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
                    master_detfoo = l1c[detector_footprint_band_name].data
                    break
            for band in self.bands:
                if band == master_band:
                    continue
                band_resolution = self.resolutions[band]
                detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
                if resolution > band_resolution:
                    factor = resolution // band_resolution
                    master_detfoo = Downsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        master_detfoo,
                        mode="masterdetfoo",
                        factor=factor,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )
                elif resolution < band_resolution:
                    factor = band_resolution // resolution
                    band_chunksize = chunksize_in_meters // band_resolution
                    resampled_detector = Upsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        mode="nearest",
                        factor=(factor, factor),
                        src_image_shape=l1c[detector_footprint_band_name].data.shape,
                        src_image_chunksize=(band_chunksize, band_chunksize),
                        depth=0,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(band_chunksize * factor, band_chunksize * factor)
                    )
                    master_detfoo = da.where(resampled_detector == master_detfoo, master_detfoo, 0)
                else:
                    resampled_detector = l1c[detector_footprint_band_name].data
                    master_detfoo = da.where(resampled_detector == master_detfoo, master_detfoo, 0)
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
                if resolution > band_resolution:
                    factor = resolution // band_resolution
                    resampled_detector = Downsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        mode="majority",
                        factor=factor,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )

                elif resolution < band_resolution:
                    factor = band_resolution // resolution
                    band_chunksize = chunksize_in_meters // band_resolution
                    resampled_detector = Upsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        mode="nearest",
                        factor=(factor, factor),
                        src_image_shape=l1c[detector_footprint_band_name].data.shape,
                        src_image_chunksize=(band_chunksize, band_chunksize),
                        depth=0,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(band_chunksize * factor, band_chunksize * factor)
                    )
                else:
                    resampled_detector = l1c[detector_footprint_band_name].data
                target_band_name = f"/conditions/mask/detector_footprint/resampled/{band}"
                resampled[target_band_name] = xr.DataArray(
                    resampled_detector,
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
            dims: Dict[str,int],
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            chunksize_in_meters: int,
            with_detfoo_filter: bool,
    ) -> xr.DataArray:
        """Adds resampled reflectance bands, returns one input band with target resolution as dummy"""
        input_band_with_target_resolution = None
        for band in self.bands:
            band_resolution = self.resolutions[band]
            source_band = l1c[f"measurements/reflectance/r{band_resolution}m/{band}"]
            factor = resolution // band_resolution
            if resolution > band_resolution:
                if downsampling == 'detectormean':
                    detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
                    resampled_band = Downsampling().apply(
                        source_band.data,
                        resampled["conditions/mask/detector_footprint/resampled/" + ("master_detfoo" if with_detfoo_filter else band)].data,
                        l1c[detector_footprint_band_name].data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=source_band.dtype,
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )
                else:
                    resampled_band = Downsampling().apply(
                        source_band.data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=source_band.dtype,
                        chunks=(chunksize_in_meters // resolution,
                                chunksize_in_meters // resolution)
                    )
            elif resolution < self.resolutions[band]:
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_band = Upsampling().apply(
                    source_band.data,
                    mode=upsampling,
                    factor=(factor, factor),
                    src_image_shape=source_band.data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    is_reflectance=True,
                    depth=overlap_depth,
                    trim=False,
                    dtype=source_band.dtype,
                    chunks=(band_chunksize * factor, band_chunksize * factor)
                )
            else:
                resampled_band = source_band.data
                if input_band_with_target_resolution is None:
                    input_band_with_target_resolution = source_band

            attrs = source_band.attrs
            attrs["coordinates"] = "crs y x lat lon"
            resampled[f"measurements/reflectance/resampled/{band}"] = xr.DataArray(
                resampled_band,
                dims=dims,
                attrs=attrs
            )
        return input_band_with_target_resolution

    def _resample_ancillary(
            self,
            ancillary: Iterable[str],
            l1c: xr.DataTree,
            lat_data: da.Array,
            lon_data: da.Array,
            dims: Dict[str, int],
            resampled: xr.DataTree
    ):
        # TODO how are ancillary bands specified?
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
            resampled[f"conditions/meteorology/emwf/{anc_band.name}_interpolated"] = xr.DataArray(anc_data, dims=dims, attrs=attrs)

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
            flag_band = f"quality/mask/r{band_resolution}m/{band}"
            if resolution > band_resolution:
                factor = resolution // band_resolution
                resampled_flags = Downsampling().apply(
                    l1c[flag_band].data,
                    mode=flagdownsampling,
                    factor=factor,
                    dtype=l1c[flag_band].dtype,
                    chunks=(chunksize_in_meters // resolution,
                            chunksize_in_meters // resolution)
                )
            elif resolution < band_resolution:
                factor = band_resolution // resolution
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_flags = Upsampling().apply(
                    l1c[flag_band].data,
                    mode='nearest',
                    factor=(factor, factor),
                    src_image_shape=l1c[flag_band].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=l1c[flag_band].dtype,
                    chunks=(band_chunksize * factor, band_chunksize * factor)
                )
            else:
                resampled_flags = l1c[flag_band].data
            attrs = l1c[flag_band].attrs.copy()
            attrs["coordinates"] = "crs y x"
            resampled[f"quality/mask/resampled/{band}"] = xr.DataArray(
                resampled_flags,
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
        band = "conditions/mask/l1c_classification/r60m/b00"
        if resolution < 60:
            factor = 60 // resolution
            band_chunksize = chunksize_in_meters // 60
            resampled_flags = Upsampling().apply(
                l1c[band].data,
                mode='nearest',
                factor=(factor, factor),
                src_image_shape=l1c[band].data.shape,
                src_image_chunksize=(band_chunksize, band_chunksize),
                depth=0,
                dtype=l1c[band].dtype,
                chunks=(band_chunksize * factor, band_chunksize * factor)
            )
        else:
            resampled_flags = l1c[band].data
            attrs = l1c[band].attrs.copy()
            attrs["coordinates"] = "crs y x"
            resampled[f"conditions/mask/l1c_classification/resampled/b00"] = xr.DataArray(
                resampled_flags,
                dims=dims,
                attrs=attrs
            )

    def _resample_sun_angles(
            self,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            input_band_with_target_resolution: xr.DataArray,
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
                input_band_with_target_resolution.data,
                tp_data=angle_data,
                resolution=resolution,
                tp_resolution=5000,
                image_shape=input_band_with_target_resolution.data.shape,
                image_chunksize=(band_chunksize, band_chunksize),
                dtype=angle_data.dtype,
            )
            resampled[f"conditions/geometry/{angle_short_name}"] = xr.DataArray(
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
        view_angle_dims = view_angle_band.dims.copy()
        view_angle_dims["x"] = dims["x"]
        view_angle_dims["y"] = dims["y"]

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
        angles_axis: int = view_angle_band.dims.index("x")
        target_chunksize = chunksize_in_meters // resolution
        detectors = view_angle_band.coords["detector"].values
        for band in self.bands:
            detector_footprint_band_name = "conditions/mask/detector_footprint/resampled/" + ("master_detfoo" if with_detfoo_filter else band)
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
            resampled_angles_data = da.stack([resampled_vza, resampled_vaa], axis=angles_axis)
            # TODO "resampled" in name?
            resampled[view_angle_band_name] = xr.DataArray(
                resampled_angles_data,
                dims=view_angle_dims,
                attrs=view_angle_band.attrs.copy(),
            )
            vza_accu.append(resampled_vza)
            vaa_accu.append(resampled_vaa)

        vza_mean = MeanAngles().apply(*vza_accu, dtype=np.float32)
        vaa_mean = MeanAngles().apply(*vaa_accu, is_azimuth_angle=True, dtype=np.float32)

        angle_mean = da.stack([vza_mean, vaa_mean], axis=1)
        resampled["conditions/geometry/mean_viewing_incidence_angles"] = xr.DataArray(
            angle_mean,
            dims={
                "band": 13,
                "angle": 2,
            }, # TODO original passed "dims", which one is correct?
        )

    #@dask.delayed
    def construct_coordinate_arrays(self, y, x, chunks: Tuple[int, int]):
        xx = da.broadcast_to(x, (y.shape[0], x.shape[0]), chunks=chunks)
        yy = da.transpose(da.broadcast_to(y, (y.shape[0], x.shape[0]), chunks=chunks[::-1]))
        return yy, xx
