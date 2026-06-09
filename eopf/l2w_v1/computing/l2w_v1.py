# Copyright 2026 ESA
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Any, Optional

import xarray as xr
import dask.array as da
import numpy as np
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from l2w_v1.exceptions.errors import MyError


class L2WV1Processor(EOProcessingUnit):
    """Polymer atmospheric correction processor.

    Wraps the hygeos/polymer run_polymer_dataset function as an
    EOPF-compliant EOProcessingUnit.
    """

    PROCESSOR_NAME = "l2w"
    PROCESSOR_VERSION = "0.0.1"
    PROCESSOR_LEVEL = "Level-2 Product"
    PROCESSOR_MODEL = False

    def __init__(self, identifier: Any = "") -> None:
        super().__init__(identifier)
        self._update_history = True

    def run(
        self,
        inputs: MappingDataType,
        adfs: Optional[MappingAuxiliary] = None,
        mode: Optional[str] = None,
        **kwargs: Any,
    ) -> MappingDataType:

        if len(inputs) < 1:
            raise MyError("Missing mandatory input product 'l1c'")

        l1c = inputs["l1c"]
        l2w_prototype = xr.open_dataset(adfs["l2w_prototype"].path, chunks=610, mask_and_scale=False) if "l2w_prototype" in adfs else None

        # introduce new L2W output
        l2w = xr.DataTree(name="l2w")

        # introduce group
        water_leaving_reflectance = xr.DataTree(name="water_leaving_reflectance")
        l2w["measurements/water_leaving_reflectance"] = water_leaving_reflectance
        shape = l1c["measurements/reflectance/r60m/b01"].data.shape

        # introduce variables
        for wavelength in [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]:
            # float Rw665(y, x);
            #         Rw665:_FillValue = NaNf;
            #         Rw665:long_name = "Water leaving reflectance at 665 nm";
            #         Rw665: units = "1";
            #         Rw665: radiation_wavelength = 665.;
            #         Rw665: wavelength_unit = "nm";
            #         Rw665: coordinates = "crs";
            band_name = f"Rw{wavelength}"
            data = l2w_prototype[band_name] if l2w_prototype else da.random.random(size=shape)
            band = xr.DataArray(
                data=data,
                dims=("y", "x"),
                attrs={
                    "long_name": f"Water leaving reflectance at {wavelength} nm",
                    "units": "1",
                    "_FillValue": np.nan,
                    "radiation_wavelength": float(wavelength),
                    "wavelength_unit": "nm",
                    "grid_mapping": "measurements/water_leaving_reflectance/crs",
                }
            )
            l2w["measurements/water_leaving_reflectance/" + band_name] = band

        # int64 crs ;
        #         string crs:crs_wkt = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         crs:semi_major_axis = 6378137. ;
        #         crs:semi_minor_axis = 6356752.31424518 ;
        #         crs:inverse_flattening = 298.257223563 ;
        #         string crs:reference_ellipsoid_name = "WGS 84" ;
        #         crs:longitude_of_prime_meridian = 0. ;
        #         string crs:prime_meridian_name = "Greenwich" ;
        #         string crs:geographic_crs_name = "WGS 84" ;
        #         string crs:horizontal_datum_name = "World Geodetic System 1984" ;
        #         string crs:projected_crs_name = "WGS 84 / UTM zone 32N" ;
        #         string crs:grid_mapping_name = "transverse_mercator" ;
        #         crs:latitude_of_projection_origin = 0. ;
        #         crs:longitude_of_central_meridian = 9. ;
        #         crs:false_easting = 500000. ;
        #         crs:false_northing = 0. ;
        #         crs:scale_factor_at_central_meridian = 0.9996 ;
        #         string crs:spatial_ref = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         string crs:GeoTransform = "399960.0 60.0 0.0 5700000.0 0.0 -60.0" ;
        #         string crs:wkt = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         string crs:i2m = "60.0,0.0,0.0,-60.0,399960.0,5700000.0" ;
        band_name = "crs"
        data = 0
        band = xr.DataArray(
            data=data,
            attrs={
                "crs_wkt": l2w_prototype["crs"].attrs["crs_wkt"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "spatial_ref": l2w_prototype["crs"].attrs["spatial_ref"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "GeoTransform": l2w_prototype["crs"].attrs["GeoTransform"] if l2w_prototype else "399960.0 60.0 0.0 5700000.0 0.0 -60.0",
                "wkt": l2w_prototype["crs"].attrs["wkt"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "i2m": l2w_prototype["crs"].attrs["i2m"] if l2w_prototype else "60.0,0.0,0.0,-60.0,399960.0,5700000.0"
            }
        )
        l2w["measurements/water_leaving_reflectance/" + band_name] = band

        band_name = "y"
        data = l1c["measurements/reflectance/r60m/y"]
        band = xr.DataArray(
            data=data,
            attrs={
                "long_name": "y coordinates of image at 60m in meters from up-left pixel",
                "units": "m"
            }
        )
        l2w["measurements/water_leaving_reflectance/" + band_name] = band

        band_name = "x"
        data = l1c["measurements/reflectance/r60m/x"]
        band = xr.DataArray(
            data=data,
            attrs={
                "long_name": "x coordinates of image at 60m in meters from up-left pixel",
                "units": "m"
            }
        )
        l2w["measurements/water_leaving_reflectance/" + band_name] = band

        # ubyte pixel_class(y, x) ;
        #          pixel_class:_FillValue = 0UB ;
        #          pixel_class:long_name = "Pixel classification flags" ;
        #          pixel_class:flag_meanings = "NO_DATA CLEAR_LAND_OR_VEGETATION CLEAR_OCEAN_WATER CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW AMBIGUOUS_CLOUD CLOUD AC_OUT_OF_BOUNDS" ;
        #          pixel_class:flag_values = 0UB, 1UB, 2UB, 3UB, 4UB, 5UB, 6UB, 7UB, 8UB, 9UB ;
        #          pixel_class:coordinates = "crs" ;
        band_name = "pixel_class"
        data = l2w_prototype[band_name] if l2w_prototype else da.random.randint(size=shape, low=0, high=10).astype(np.ubyte)
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Pixel classification flags",
                "flag_values": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "flag_meanings": "NO_DATA CLEAR_LAND_OR_VEGETATION CLEAR_OCEAN_WATER CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW AMBIGUOUS_CLOUD CLOUD AC_OUT_OF_BOUNDS",
                "_FillValue": np.ubyte(0),
                "grid_mapping": "quality/water/crs",
            }
        )
        l2w["quality/water/" + band_name] = band

        # ubyte sen2water_flags(y, x) ;
        #         string sen2water_flags:long_name = "quality and algorithm flags" ;
        #         string sen2water_flags:flag_meanings = "c2rcc_oor acolite_negatives polymer_invalid with_c2rcc with_acolite with_polymer" ;
        #         sen2water_flags:flag_masks = 1UB, 2UB, 4UB, 8UB, 16UB, 32UB ;
        #         string sen2water_flags:coordinates = "crs" ;
        band_name = "sen2water_flags"
        data = l2w_prototype[band_name] if l2w_prototype else da.random.randint(size=shape, low=0, high=32).astype(np.ubyte)
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Quality and algorithm flags",
                "flag_masks": [1, 2, 4, 8, 16, 32],
                "flag_meanings": "c2rcc_oor acolite_negatives polymer_invalid with_c2rcc with_acolite with_polymer",
                "grid_mapping": "quality/water/crs",
            }
        )
        l2w["quality/water/" + band_name] = band

        # int pixel_classif_flags(y, x) ;
        #          string pixel_classif_flags:flag_meanings = "IDEPIX_INVALID IDEPIX_CLOUD IDEPIX_CLOUD_AMBIGUOUS IDEPIX_CLOUD_SURE IDEPIX_CLOUD_BUFFER IDEPIX_CLOUD_SHADOW IDEPIX_SNOW_ICE IDEPIX_BRIGHT IDEPIX_WHITE IDEPIX_COASTLINE IDEPIX_LAND IDEPIX_CIRRUS_SURE IDEPIX_CIRRUS_AMBIGUOUS IDEPIX_CLEAR_LAND IDEPIX_CLEAR_WATER IDEPIX_WATER IDEPIX_BRIGHTWHITE IDEPIX_VEG_RISK IDEPIX_MOUNTAIN_SHADOW IDEPIX_POTENTIAL_SHADOW IDEPIX_CLUSTERED_CLOUD_SHADOW" ;
        #          pixel_classif_flags:flag_masks = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 ;
        #          string pixel_classif_flags:flag_coding_name = "pixel_classif_flags" ;
        #          string pixel_classif_flags:flag_descriptions = "Invalid pixels\tPixels which are either cloud_sure or cloud_ambiguous\tSemi transparent clouds, or clouds where the detection level is uncertain\tFully opaque clouds with full confidence of their detection\tA buffer of n pixels around a cloud. n is a user supplied parameter. Applied to pixels masked as \'cloud\'\tPixel is affected by a cloud shadow (combination of shifted cloud mask in cloud gaps and dark clusters coinciding with a corrected shifted cloud mask)\tClear snow/ice pixels\tBright pixels\tWhite pixels\tPixels at a coastline\tLand pixels\tCirrus clouds with full confidence of their detection\tCirrus clouds, or clouds where the detection level is uncertain\tClear land pixels\tClear water pixels\tWater pixels\t\'Brightwhite\' pixels\tPixels with vegetation risk\tPixel is affected by mountain shadow\tPotentially a cloud shadow pixel\tCloud shadow identified by clustering algorithm" ;
        #          pixel_class:coordinates = "crs" ;
        band_name = "pixel_classif_flags"
        data = l2w_prototype[band_name] if l2w_prototype else da.random.randint(size=shape, low=0, high=2*1048576)
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Idepix pixel classification flags",
                "flag_values": [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576],
                "flag_meanings": "IDEPIX_INVALID IDEPIX_CLOUD IDEPIX_CLOUD_AMBIGUOUS IDEPIX_CLOUD_SURE IDEPIX_CLOUD_BUFFER IDEPIX_CLOUD_SHADOW IDEPIX_SNOW_ICE IDEPIX_BRIGHT IDEPIX_WHITE IDEPIX_COASTLINE IDEPIX_LAND IDEPIX_CIRRUS_SURE IDEPIX_CIRRUS_AMBIGUOUS IDEPIX_CLEAR_LAND IDEPIX_CLEAR_WATER IDEPIX_WATER IDEPIX_BRIGHTWHITE IDEPIX_VEG_RISK IDEPIX_MOUNTAIN_SHADOW IDEPIX_POTENTIAL_SHADOW IDEPIX_CLUSTERED_CLOUD_SHADOW",
                "_FillValue": 0,
                "grid_mapping": "quality/water/crs",
            }
        )
        l2w["quality/water/" + band_name] = band

        # int64 crs ;
        #         string crs:crs_wkt = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         crs:semi_major_axis = 6378137. ;
        #         crs:semi_minor_axis = 6356752.31424518 ;
        #         crs:inverse_flattening = 298.257223563 ;
        #         string crs:reference_ellipsoid_name = "WGS 84" ;
        #         crs:longitude_of_prime_meridian = 0. ;
        #         string crs:prime_meridian_name = "Greenwich" ;
        #         string crs:geographic_crs_name = "WGS 84" ;
        #         string crs:horizontal_datum_name = "World Geodetic System 1984" ;
        #         string crs:projected_crs_name = "WGS 84 / UTM zone 32N" ;
        #         string crs:grid_mapping_name = "transverse_mercator" ;
        #         crs:latitude_of_projection_origin = 0. ;
        #         crs:longitude_of_central_meridian = 9. ;
        #         crs:false_easting = 500000. ;
        #         crs:false_northing = 0. ;
        #         crs:scale_factor_at_central_meridian = 0.9996 ;
        #         string crs:spatial_ref = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         string crs:GeoTransform = "399960.0 60.0 0.0 5700000.0 0.0 -60.0" ;
        #         string crs:wkt = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]" ;
        #         string crs:i2m = "60.0,0.0,0.0,-60.0,399960.0,5700000.0" ;
        band_name = "crs"
        data = 0
        band = xr.DataArray(
            data=data,
            attrs={
                "crs_wkt": l2w_prototype["crs"].attrs["crs_wkt"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "spatial_ref": l2w_prototype["crs"].attrs["spatial_ref"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "GeoTransform": l2w_prototype["crs"].attrs["GeoTransform"] if l2w_prototype else "399960.0 60.0 0.0 5700000.0 0.0 -60.0",
                "wkt": l2w_prototype["crs"].attrs["wkt"] if l2w_prototype else "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]",
                "i2m": l2w_prototype["crs"].attrs["i2m"] if l2w_prototype else "60.0,0.0,0.0,-60.0,399960.0,5700000.0"
            }
        )
        l2w["quality/water/" + band_name] = band

        band_name = "y"
        data = l1c["measurements/reflectance/r60m/y"]
        band = xr.DataArray(
            data=data,
            attrs={
                "long_name": "y coordinates of image at 60m in meters from up-left pixel",
                "units": "m"
            }
        )
        l2w["quality/water/" + band_name] = band

        band_name = "x"
        data = l1c["measurements/reflectance/r60m/x"]
        band = xr.DataArray(
            data=data,
            attrs={
                "long_name": "x coordinates of image at 60m in meters from up-left pixel",
                "units": "m"
            }
        )
        l2w["quality/water/" + band_name] = band

        l2w.attrs["stac_discovery"] = {
            "type": "Feature",
            "id": "S02MSIL2W_20210627T100559_0000_B022_TFBD",
            "bbox": [
                11.755584916064224,
                45.93349650488509,
                10.290484068727272,
                46.94616825397537,
            ],
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [10.314035280374382, 46.94616825397537],
                        [11.755584916064224, 46.920536388192055],
                        [11.706254333045143, 45.93349650488509],
                        [10.290484068727272, 45.95826451839171],
                        [10.314035280374382, 46.94616825397537],
                    ]
                ],
            },
            "properties": {
                "product:type": "S02MSIL2W",
                "processing:facility": "ESA",
                "processing:level": "L2W",
                "processing:lineage": "systematic",
                "processing:version": "00.08",
                "product:timeliness": "PT3H",
                "product:timeliness_category": "NRT",
                "constellation": "sentinel-2",
                "platform": "sentinel-2b",
                "sat:platform_international_designator": "2015-028A",
                "instruments": ["msi"],
                "sat:absolute_orbit": 22499,
                "sat:orbit_state": "descending",
                "sat:relative_orbit": 22,
                "eopf:datastrip_id": "S2B_OPER_MSI_L1C_DS_S2RP_20230318T143834_S20210627T101702_N05.00",
                "eopf:datatake_id": "GS2B_20210627T100559_022499_N05.00",
                "eopf:instrument_mode": "INS-NOBS",
                "gsd": "60",
                "proj:bbox": [600000.0, 5090220.0, 709800.0, 5200020.0],
                "proj:code": "EPSG:32632",
                "proj:transform": [
                    60.0,
                    0.0,
                    600000.0,
                    0.0,
                    -60.0,
                    5200020.0,
                    0.0,
                    0.0,
                    1.0,
                ],
                "proj:wkt2": 'PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32632"]]',
                "start_datetime": "2021-06-27T10:05:59.024Z",
                "end_datetime": "2021-06-27T10:05:59.024Z",
                "datetime": "2021-06-27T10:05:59.024000Z",
                "eo:cloud_cover": 15.3870823255397,
                "eo:snow_cover": 8.31696975126161,
                "providers": [{"name": "ESA", "roles": ["producer"]}],
                "created": "2023-03-18T14:38:34.000000Z",
                "sci:doi": "TBD",
                "bands": [
                    {
                        "eo:center_wavelength": 0.4423,
                        "eo:full_width_half_max": 0.02,
                        "name": "Rw443",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.4923,
                        "eo:full_width_half_max": 0.065,
                        "name": "Rw490",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.559,
                        "eo:full_width_half_max": 0.035,
                        "name": "Rw560",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.665,
                        "eo:full_width_half_max": 0.03,
                        "name": "Rw665",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.7038,
                        "eo:full_width_half_max": 0.015,
                        "name": "Rw703",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.7391000000000001,
                        "eo:full_width_half_max": 0.015,
                        "name": "Rw740",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.7797000000000001,
                        "eo:full_width_half_max": 0.02,
                        "name": "Rw783",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.833,
                        "eo:full_width_half_max": 0.105,
                        "name": "Rw842",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.864,
                        "eo:full_width_half_max": 0.02,
                        "name": "Rw865",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 0.9432,
                        "eo:full_width_half_max": 0.02,
                        "name": "Rw945",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 1.3769,
                        "eo:full_width_half_max": 0.03,
                        "name": "Rw1375",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 1.6104,
                        "eo:full_width_half_max": 0.09,
                        "name": "Rw1610",
                        "raster:spatial_resolution": "60",
                    },
                    {
                        "eo:center_wavelength": 2.1856999999999998,
                        "eo:full_width_half_max": 0.18,
                        "name": "Rw2190",
                        "raster:spatial_resolution": "60",
                    },
                ],
            },
            "assets": {
                "Rw443": {
                    "title": "Rw443",
                    "href": "/measurements/water_leaving_reflectance/Rw443",
                },
                "Rw490": {
                    "title": "Rw490",
                    "href": "/measurements/water_leaving_reflectance/Rw490",
                },
                "Rw560": {
                    "title": "Rw560",
                    "href": "/measurements/water_leaving_reflectance/Rw560",
                },
                "Rw665": {
                    "title": "Rw665",
                    "href": "/measurements/water_leaving_reflectance/Rw665",
                },
                "Rw705": {
                    "title": "Rw705",
                    "href": "/measurements/water_leaving_reflectance/Rw705",
                },
                "Rw740": {
                    "title": "Rw740",
                    "href": "/measurements/water_leaving_reflectance/Rw740",
                },
                "Rw783": {
                    "title": "Rw783",
                    "href": "/measurements/water_leaving_reflectance/Rw783",
                },
                "Rw842": {
                    "title": "Rw842",
                    "href": "/measurements/water_leaving_reflectance/Rw842",
                },
                "Rw865": {
                    "title": "Rw865",
                    "href": "/measurements/water_leaving_reflectance/Rw865",
                },
                "Rw945": {
                    "title": "Rw945",
                    "href": "/measurements/water_leaving_reflectance/Rw945",
                },
                "Rw1375": {
                    "title": "Rw1375",
                    "href": "/measurements/water_leaving_reflectance/Rw1375",
                },
                "Rw1610": {
                    "title": "Rw1610",
                    "href": "/measurements/water_leaving_reflectance/Rw1610",
                },
                "Rw2190": {
                    "title": "Rw2190",
                    "href": "/measurements/water_leaving_reflectance/Rw2190",
                },
                "pixel_class": {
                    "title": "pixel_class",
                    "href": "/quality/water/pixel_class",
                },
                "sen2water_flags": {
                    "title": "sen2water_flags",
                    "href": "/quality/water/sen2water_flags",
                },
                "pixel_classif_flags": {
                    "title": "pixel_classif_flags",
                    "href": "/quality/water/pixel_classif_flags",
                },
            },
            "stac_version": "1.1.0",
            "stac_extensions": [
                "https://stac-extensions.github.io/eopf/v1.2.0/schema.json",
                "https://stac-extensions.github.io/eo/v1.1.0/schema.json",
                "https://stac-extensions.github.io/raster/v2.0.0/schema.json",
                "https://stac-extensions.github.io/sat/v1.1.0/schema.json",
                "https://stac-extensions.github.io/view/v1.0.0/schema.json",
                "https://stac-extensions.github.io/scientific/v1.0.0/schema.json",
                "https://stac-extensions.github.io/processing/v1.2.0/schema.json",
                "https://stac-extensions.github.io/product/v0.1.0/schema.json",
            ],
        }
        l2w.attrs["other_metadata" ] = {
            "title": "OPT-MPC Sen2Water water-leaving reflectances in 60m",
            "tracking_id": "21ccb604-ba3d-11f0-b1d5-0433c22b8076",
            "pixel_statistics": {
                "ocean": {
                    "valid_count": 213224,
                    "clear_count": 140090,
                    "snow_ice_count": 2,
                    "cloud_count": 73132,
                },
                "inland_water": {
                    "valid_count": 1783,
                    "clear_count": 1031,
                    "snow_ice_count": 0,
                    "cloud_count": 752,
                },
                "land": {
                    "valid_count": 225938,
                    "clear_count": 161244,
                    "snow_ice_count": 3,
                    "cloud_count": 64691,
                },
                "valid_count": 440945,
            },
            "eopf_category": "eoproduct",
        }

        outputs = {"l2w": l2w}

        l2w.attrs["processing_history"] = l1c.attrs["processing_history"]
        # TODO seems to happen outside of the processor, lacks DEM and S2W water mask
        # self._add_history_event(inputs, outputs=outputs, mode="default", adfs=adfs)

        # return inputs
        return outputs
