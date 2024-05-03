# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

from typing import Tuple, Union, Dict, List

import dask.array as da
import xarray as xr
import numpy as np
from pyproj import Transformer, CRS
from sen2water.eoutils.eoprocessingifc import Operator
from sen2water.eoutils.eoutils import copy_variable
from sen2water.msiresampling.ancillaryinterpolation import AncillaryInterpolation
from sen2water.msiresampling.downsampling import Downsampling
from sen2water.msiresampling.geocoordinates import GeoCoordinates
from sen2water.msiresampling.meanangles import MeanAngles
from sen2water.msiresampling.upsampling import Upsampling
from sen2water.msiresampling.tpinterpolation import TpInterpolation
from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.constants import MsiConstants


class ResamplingOperator(Operator):

    bands = MsiConstants.bands
    resolutions = MsiConstants.resolutions
    flag_band_prefixes = MsiConstants.flag_band_prefixes
    cloud_ice_flags = MsiConstants.cloud_ice_flags
    ancillary_bands = MsiConstants.ancillary_bands
    sun_angles = MsiConstants.sun_angles
    detector_colour = MsiConstants.detector_colour
    quality_flags_colour = MsiConstants.quality_flags_colour
    cloud_ice_colour = MsiConstants.cloud_ice_colour

    chunksize_in_meters = 36600

    def preferred_chunks(self):
        return { "y10": self.chunksize_in_meters // 10,
                 "x10": self.chunksize_in_meters // 10,
                 "y20": self.chunksize_in_meters // 20,
                 "x20": self.chunksize_in_meters // 20,
                 "y60": self.chunksize_in_meters // 60,
                 "x60": self.chunksize_in_meters // 60,
                 }

    def run(
        self,
        l1c: xr.Dataset,
        resolution: int,
        downsampling: Union['detectormean', 'first', 'min', 'max', 'mean', 'median'],
        flagdownsampling: str,
        upsampling: Union["nearest", "bilinear", "bicubic"],
        ancillary: List[str],
        with_master_detfoo: bool,
        merge_flags: bool=False
    ) -> xr.Dataset:

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {"y": l1c.sizes[f"y{resolution}"], "x": l1c.sizes[f"x{resolution}"]}
        chunks = (self.chunksize_in_meters // resolution, self.chunksize_in_meters // resolution)
        coordinate_bands = {}
        target_bands = {}

        # copy 1-d metric coordinate bands

        coordinate_bands["y"] = xr.DataArray(
            l1c[f"y{resolution}"].data,
            dims=['y'],
            attrs={"long_name": "y coordinate of projection",
                   "standard_name": "projection_y_coordinate",
                   "units": "m"}
        )
        coordinate_bands["x"] = xr.DataArray(
            l1c[f"x{resolution}"].data,
            dims=['x'],
            attrs={"long_name": "x coordinate of projection",
                   "standard_name": "projection_x_coordinate",
                   "units": "m"}
        )

        # transform UTM coordinates into geo-coordinates, add them

        y = da.from_array(l1c[f"y{resolution}"].values, chunks=self.chunksize_in_meters // resolution)
        x = da.from_array(l1c[f"x{resolution}"].values, chunks=self.chunksize_in_meters // resolution)
        xx = da.broadcast_to(x, (y.shape[0], x.shape[0]), chunks=chunks)
        yy = da.transpose(da.broadcast_to(y, (y.shape[0], x.shape[0]), chunks=chunks[::-1]))
        transformer = Transformer.from_crs(CRS.from_cf(l1c['spatial_ref'].attrs), CRS("EPSG:4326"))
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

        coordinate_bands["lat"] = xr.DataArray(
            lat_data,
            dims=['y', 'x'],
            attrs={"standard_name": "latitude",
                   "units": "degrees_north"}
        )
        coordinate_bands["lon"] = xr.DataArray(
            lon_data,
            dims=['y', 'x'],
            attrs={"standard_name": "longitude",
                   "units": "degrees_east"}
        )

        coordinate_bands["crs"] = copy_variable(l1c[f"spatial_ref_{resolution}"])
        # make SNAP happy
        if "wkt" not in coordinate_bands["crs"].attrs:
            coordinate_bands["crs"].attrs["wkt"] = coordinate_bands["crs"].attrs["crs_wkt"]
        # GeoTransform 699960.0 60.0 0.0 5200020.0 0.0 -60.0
        # i2m 60.0,0.0,0.0,-60.0,699960.0,5200020.0
        if "i2m" not in coordinate_bands["crs"].attrs:
            geo_t = coordinate_bands["crs"].attrs["GeoTransform"].split(" ")
            coordinate_bands["crs"].attrs["i2m"] = ",".join([geo_t[i] for i in [1,2,4,5,0,3]])

        # copy geographic coordinates of ancillary data

        for band in ["aux_latitude", "aux_longitude"]:
            attrs = dict(l1c[band].attrs)
            attrs["offset_y"] = 0.0
            attrs["offset_x"] = 0.0
            attrs["subsampling_y"] = 10980 * 10 / resolution / 8
            attrs["subsampling_x"] = 10980 * 10 / resolution / 8
            attrs["_FillValue"] = np.nan
            coordinate_bands[band] = xr.DataArray(l1c[band].data,
                                                  dims={f"tp_{band[4:7]}": l1c[band].sizes[band]},
                                                  attrs=attrs)

        # select detector per target pixel (with_master_detfoo) or per target pixel per band

        self._resample_detectors(resolution, dims, l1c, target_bands, with_master_detfoo=with_master_detfoo)

        # resample reflectance bands B1 .. B12

        input_band_with_target_resolution = self._resample_reflectance(
            downsampling, upsampling, resolution, overlap_depth, dims, l1c, target_bands, with_master_detfoo=with_master_detfoo)

        # resample flag bands B_xxx_B1 .. B_zzz_B12

        self._resample_flags(flagdownsampling, resolution, overlap_depth, dims, l1c, target_bands, merge_flags)

        # resample cloud and ice flag bands B_opaque_clouds, ...

        self._resample_cloud_ice(flagdownsampling, resolution, dims, l1c, target_bands, merge_flags)

        # interpolate ancillary bands if requested

        if ancillary:
            self._resample_ancillary(ancillary, l1c, lat_data, lon_data, dims, target_bands)
        del lat_data, lon_data

        # copy ECMWF and CAMS bands tco3, ...

        for band in self.ancillary_bands:
            attrs = {n: l1c[band].attrs[n] for n in l1c[band].attrs if not n.startswith("GRIB_") and n != "_FillValue"}
            attrs["coordinates"] = "tp_lat tp_lon"
            attrs["offset_y"] = 0.0
            attrs["offset_x"] = 0.0
            attrs["subsampling_y"] = 10980 * 10 / resolution / 8
            attrs["subsampling_x"] = 10980 * 10 / resolution / 8
            tp_dims = {"tp_lat": l1c[band].sizes["aux_latitude"],
                       "tp_lon": l1c[band].sizes["aux_longitude"]}
            target_bands[band] = xr.DataArray(l1c[band].data,
                                              dims=tp_dims,
                                              attrs=attrs)

        # interpolate sun angles

        self._resample_sun_angles(resolution, dims, input_band_with_target_resolution, l1c, target_bands)

        # interpolate viewing angles considering detector footprint

        self._resample_viewing_angles(resolution, dims, l1c, target_bands, with_master_detfoo=with_master_detfoo)

        self._add_snap_masks(merge_flags, target_bands, with_master_detfoo=with_master_detfoo)

        # for debugging copy angles per band and detector
        # for band in self.bands:
        #     for detector_i, detector in enumerate(detectors):
        #         target_bands[f"vaa_{band}_{detector}"] = xr.DataArray(
        #             np.copy(l1c[f"vaa_{band}"].values[detector_i]),
        #             dims={"angle_y": 23, "angle_x": 23},
        #             attrs=l1c[f"vaa_{band}"].attrs)

        # create attributes similar to S2Resampling
        # add spacecraft to global attributes because we do not provide the metadata variable

        target_attrs = {
            "conventions": "CF-1.10",
            "TileSize": "610:610",
            "product_type": "S2_MSI_Level-1C",
            "platform": l1c.attrs["spacecraft"],
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
            "start_date": l1c.attrs["start_date"],
            "stop_date": l1c.attrs["stop_date"],
        };

        return xr.Dataset(target_bands, coordinate_bands, target_attrs)

    def _resample_ancillary(self, ancillary, l1c, lat_data, lon_data, dims, target_bands):
        for anc_band in ancillary:
            anc_data = AncillaryInterpolation().apply(
                lat_data,
                lon_data,
                anc_lat=l1c["aux_latitude"].values,
                anc_lon=l1c["aux_longitude"].values,
                anc_data=l1c[anc_band].values,
                variable=anc_band,
                dtype=np.float32,
            )
            attrs = {n: l1c[anc_band].attrs[n] for n in l1c[anc_band].attrs if
                     not n.startswith("GRIB_") and n != "_FillValue"}
            attrs["coordinates"] = "crs y x lat lon"
            target_bands[f"{anc_band}_interpolated"] = xr.DataArray(anc_data, dims=dims, attrs=attrs)

    def _resample_reflectance(
            self,
            downsampling: Union[
                "detectormean", "mean", "median", "min", "max", "first",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
            ],
            upsampling: Union["nearest", "bilinear", "bicubic"],
            resolution: int,
            overlap_depth: int,
            dims: Dict[str,int],
            l1c: xr.Dataset,
            target_bands: Dict[str, xr.DataArray],
            with_master_detfoo: bool
    ) -> xr.DataArray:
        """Adds resampled reflectance bands, returns one input band with target resolution as dummy"""
        input_band_with_target_resolution = None
        for band in self.bands:
            if resolution > self.resolutions[band]:
                factor = resolution // self.resolutions[band]
                if downsampling == 'detectormean':
                    detector_footprint_band_name = f"B_detector_footprint_{band}"
                    resampled = Downsampling().apply(
                        l1c[band].data,
                        target_bands["master_detfoo" if with_master_detfoo else detector_footprint_band_name].data,
                        l1c[detector_footprint_band_name].data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=l1c[band].dtype,
                        chunks=(self.chunksize_in_meters // resolution,
                                self.chunksize_in_meters // resolution)
                    )
                else:
                    resampled = Downsampling().apply(
                        l1c[band].data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=l1c[band].dtype,
                        chunks=(self.chunksize_in_meters // resolution,
                                self.chunksize_in_meters // resolution)
                    )
            elif resolution < self.resolutions[band]:
                factor = self.resolutions[band] // resolution
                band_chunksize = self.chunksize_in_meters // self.resolutions[band]
                resampled = Upsampling().apply(
                    l1c[band].data,
                    mode=upsampling,
                    factor=(factor, factor),
                    src_image_shape=l1c[band].data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    is_reflectance=True,
                    depth=overlap_depth,
                    trim=False,
                    dtype=l1c[band].dtype,
                    chunks=(band_chunksize * factor, band_chunksize * factor)
                )
            else:
                resampled = l1c[band].data
                if input_band_with_target_resolution is None:
                    input_band_with_target_resolution = l1c[band]
            attrs = l1c[band].attrs
            attrs["coordinates"] = "crs y x lat lon"
            target_bands[band] = xr.DataArray(
                resampled,
                dims=dims,
                attrs=attrs
            )
        return input_band_with_target_resolution

    def _resample_flags(
            self,
            flagdownsampling: Union[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            overlap_depth: int,
            dims: Dict[str, int],
            l1c: xr.Dataset,
            target_bands: Dict[str, xr.DataArray],
            merge_flags: bool
    ):
        """Adds quality flags bands in target resolution"""
        for band in self.bands:
            for flag_band_prefix in (["quality_flags"] if merge_flags else self.flag_band_prefixes):
                flag_band = f"{flag_band_prefix}_{band}"
                if resolution > self.resolutions[band]:
                    factor = resolution // self.resolutions[band]
                    resampled = Downsampling().apply(
                        l1c[flag_band].data,
                        mode=flagdownsampling,
                        factor=factor,
                        dtype=l1c[flag_band].dtype,
                        chunks=(self.chunksize_in_meters // resolution,
                                self.chunksize_in_meters // resolution)
                    )
                elif resolution < self.resolutions[band]:
                    factor = self.resolutions[band] // resolution
                    band_chunksize = self.chunksize_in_meters // self.resolutions[band]
                    resampled = Upsampling().apply(
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
                    resampled = l1c[flag_band].data
                attrs = {"long_name": f"quality mask of {band}",
                         "coordinates": "crs y x lat lon"}
                if merge_flags:
                    attrs["flag_masks"] = l1c[flag_band].attrs["flag_masks"]
                    attrs["flag_meanings"] = l1c[flag_band].attrs["flag_meanings"]
                target_bands[flag_band] = xr.DataArray(
                    resampled,
                    dims=dims,
                    attrs=attrs
                )

    def _resample_cloud_ice(
            self,
            flagdownsampling: Union[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            dims: Dict[str, int],
            l1c: xr.Dataset,
            target_bands: Dict[str, xr.DataArray],
            merge_flags: bool
    ):
        """Adds cloud and ice flags in target resolution"""
        for band in (["cloud_ice_flags"] if merge_flags else self.cloud_ice_flags):
            if resolution < 60:
                factor = 60 // resolution
                band_chunksize = self.chunksize_in_meters // 60
                resampled = Upsampling().apply(
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
                resampled = l1c[band].data
            attrs = {"long_name": f"quality classification",
                     "coordinates": "crs y x lat lon"}
            if merge_flags:
                attrs["flag_masks"] = l1c[band].attrs["flag_masks"]
                attrs["flag_meanings"] = l1c[band].attrs["flag_meanings"]
            target_bands[band] = xr.DataArray(
                resampled,
                dims=dims,
                attrs=attrs
            )

    def _resample_sun_angles(
            self,
            resolution: int,
            dims: Tuple[int, int],
            input_band_with_target_resolution: xr.DataArray,
            l1c: xr.Dataset,
            target_bands: Dict[str,xr.DataArray]
    ):
        """Adds sun angles resampled from TP grid in target resolution"""
        for band in self.sun_angles:
            band_chunksize = self.chunksize_in_meters // resolution
            resampled = TpInterpolation().apply(
                input_band_with_target_resolution.data,
                tp_data=l1c[band].data,
                resolution=resolution,
                tp_resolution=5000,
                image_shape=input_band_with_target_resolution.data.shape,
                image_chunksize=(band_chunksize, band_chunksize),
                dtype=l1c[band].dtype,
            )
            target_bands[band] = xr.DataArray(
                resampled,
                dims=dims,
                attrs={"long_name": f"Solar {band[4:]} angle",
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x lat lon"}
            )

    def _resample_detectors(
            self,
            resolution: int,
            dims: Dict[str,int],
            l1c: xr.Dataset,
            target_bands: Dict[str,xr.DataArray],
            with_master_detfoo: bool=True
    ):
        """Adds resampled detector index"""
        if with_master_detfoo:
            # The result shall be a single detfoo with fill value in the overlapping area without a common detfoo value.
            # Any of the source bands with the target resolution determines the master detector footprint value.
            # If any of the other bands does not have a contribution with this master detfoo value then the master detfoo pixel is set to invalid.
            master_band = None
            master_detfoo = None
            for band in self.bands:
                if self.resolutions[band] == resolution:
                    master_band = band
                    detector_footprint_band_name = f"B_detector_footprint_{band}"
                    master_detfoo = l1c[detector_footprint_band_name].data
                    break
            for band in self.bands:
                if band == master_band:
                    continue
                detector_footprint_band_name = f"B_detector_footprint_{band}"
                band_resolution = self.resolutions[band]
                if resolution > band_resolution:
                    factor = resolution // band_resolution
                    master_detfoo = Downsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        master_detfoo,
                        mode="masterdetfoo",
                        factor=factor,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(self.chunksize_in_meters // resolution,
                                self.chunksize_in_meters // resolution)
                    )
                elif resolution < band_resolution:
                    factor = band_resolution // resolution
                    band_chunksize = self.chunksize_in_meters // band_resolution
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
            target_bands["master_detfoo"] = xr.DataArray(
                master_detfoo,
                dims=dims,
                attrs={"long_name": f"detector footprint of {band}",
                       "_FillValue": np.uint8(0),
                       "coordinates": "crs y x lat lon"}
            )
        else:
            for band in self.bands:
                detector_footprint_band_name = f"B_detector_footprint_{band}"
                band_resolution = self.resolutions[band]
                if resolution > band_resolution:
                    factor = resolution // band_resolution
                    resampled_detector = Downsampling().apply(
                        l1c[detector_footprint_band_name].data,
                        mode="majority",
                        factor=factor,
                        dtype=l1c[detector_footprint_band_name].dtype,
                        chunks=(self.chunksize_in_meters // resolution,
                                self.chunksize_in_meters // resolution)
                    )

                elif resolution < band_resolution:
                    factor = band_resolution // resolution
                    band_chunksize = self.chunksize_in_meters // band_resolution
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
                target_bands[detector_footprint_band_name] = xr.DataArray(
                    resampled_detector,
                    dims=dims,
                    attrs={"long_name": f"detector footprint of {band}",
                           "_FillValue": np.uint8(0),
                           "coordinates": "crs y x lat lon"}
                )

    def _resample_viewing_angles(
            self,
            resolution: int,
            dims: Dict[str,int],
            l1c: xr.Dataset,
            target_bands: Dict[str,xr.DataArray],
            with_master_detfoo: bool=True
    ):
        """Resamples viewing angles per detector and adds viewing angles per band"""
        vza_accu = []
        vaa_accu = []
        for band in self.bands:
            detector_footprint_band_name = "master_detfoo" if with_master_detfoo else f"B_detector_footprint_{band}"
            vza_band_name = f"view_zenith_{band}"
            vaa_band_name = f"view_azimuth_{band}"
            target_chunksize = self.chunksize_in_meters // resolution
            detectors = l1c["detector"].values
            detector_footprint = target_bands[detector_footprint_band_name].data
            extended_vza, extended_vaa = AnglesInterpolation().expand_angles_per_detector(
                l1c[f"vza_{band}"].values,
                l1c[f"vaa_{band}"].values
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
            target_bands[vza_band_name] = xr.DataArray(
                resampled_vza,
                dims=dims,
                attrs={"long_name": "Viewing incidence zenith angle",
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x lat lon"}
            )
            target_bands[vaa_band_name] = xr.DataArray(
                resampled_vaa,
                dims=dims,
                attrs={"long_name": "Viewing incidence azimuth angle",
                       "units": "degrees",
                       "_FillValue": np.nan,
                       "coordinates": "crs y x lat lon"}
            )
            vza_accu.append(resampled_vza)
            vaa_accu.append(resampled_vaa)

        vza_mean = MeanAngles().apply(*vza_accu, dtype=np.float32)
        target_bands["view_zenith_mean"] = xr.DataArray(
            vza_mean,
            dims=dims,
            attrs={"long_name": "Viewing incidence zenith angle",
                   "units": "degrees",
                   "_FillValue": np.nan,
                   "coordinates": "crs y x lat lon"}
        )
        vaa_mean = MeanAngles().apply(*vaa_accu, is_azimuth_angle=True, dtype=np.float32)
        target_bands["view_azimuth_mean"] = xr.DataArray(
            vaa_mean,
            dims=dims,
            attrs={"long_name": "Viewing incidence azimuth angle",
                   "units": "degrees",
                   "_FillValue": np.nan,
                   "coordinates": "crs y x lat lon"}
        )

    def _add_snap_masks(
            self,
            merge_flags: bool,
            target_bands: Dict[str,xr.DataArray],
            with_master_detfoo: bool=True
    ):
        """Adds mask expressions without extend for SNAP"""
        if with_master_detfoo:
            for detector in range(1, 13):
                ResamplingOperator._add_snap_mask(f"detector_footprint-{detector}_mask",
                                    f"master_detfoo=={detector}",
                                    [np.int32(c) for c in self.detector_colour[detector-1]],
                                    0.5,
                                    target_bands)
        else:
            for band in self.bands:
                for detector in range(1, 13):
                    ResamplingOperator._add_snap_mask(f"detector_footprint-{band}-{detector}_mask",
                                        f"B_detector_footprint_{band}=={detector}",
                                        [np.int32(c) for c in self.detector_colour[detector-1]],
                                        0.5,
                                        target_bands)
        for band in self.bands:
            for flag_i, flag in enumerate(self.flag_band_prefixes):
                ResamplingOperator._add_snap_mask(f"{flag[2:]}_{band}_mask",
                                    f"quality_flags_{band} & {1<<flag_i} != 0" if merge_flags else f"{flag}_{band} != 0",
                                    [np.int32(c) for c in self.quality_flags_colour[flag_i]],
                                    0.5,
                                    target_bands)
        for flag_i, flag in enumerate(self.cloud_ice_flags):
            ResamplingOperator._add_snap_mask(f"{flag}_mask",
                                f"cloud_ice_flags=={1<<flag_i} != 0" if merge_flags else f"{flag} != 0",
                                [np.int32(c) for c in self.cloud_ice_colour[flag_i]],
                                0.5,
                                target_bands)

    @staticmethod
    def _add_snap_mask(
            name: str,
            expression: str,
            colours: Tuple[int,int,int,int],
            transparency: float,
            target_bands: Dict[str,xr.DataArray]
    ):
        """Creates one mask expression as band without dims"""
        target_bands[name] = xr.DataArray(
            np.byte(0),
            attrs={
                "title": name,
                "expression": expression,
                "color": colours,
                "transparency": transparency
            }
        )
