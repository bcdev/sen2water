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
from eopf.logging import EOLogging

from l2w_v1.exceptions.errors import MyError


class L2WV1Processor(EOProcessingUnit):
    """Polymer atmospheric correction processor.

    Wraps the hygeos/polymer run_polymer_dataset function as an
    EOPF-compliant EOProcessingUnit.
    """

    PROCESSOR_NAME = "l2w"
    PROCESSOR_VERSION = "0.0.1"
    PROCESSOR_LEVEL = "L2"
    PROCESSOR_MODEL = False

    def run(
        self,
        inputs: MappingDataType,
        adfs: Optional[MappingAuxiliary] = None,
        mode: Optional[str] = None,
        **kwargs: Any,
    ) -> MappingDataType:

        logger = EOLogging().get_logger()

        if len(inputs) < 1:
            raise MyError("Missing mandatory input product 'l1c'")

        l1c = inputs["l1c"]
        l2w_prototype = adfs["l2w_prototype"] if "l2w_prototype" in adfs else None

        # introduce group
        water_leaving_reflectance = xr.DataTree(name="water_leaving_reflectance")
        l1c["measurements/water_leaving_reflectance"] = water_leaving_reflectance
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
            data = l2w_prototype[band_name] if l2w_prototype else da.random.random(size=shape, dtype=float)
            band = xr.DataArray(
                data=data,
                dims=("y", "x"),
                attrs={
                    "long_name": f"Water leaving reflectance at {wavelength} nm",
                    "units": "1",
                    "_FillValue": np.nanf,
                    "radiation_wavelength": float(wavelength),
                    "wavelength_unit": "nm",
                    "grid_mapping": "crs",
                }
            )
            l1c["measurements/water_leaving_reflectance/" + band_name] = band

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
        l1c["measurements/water_leaving_reflectance/" + band_name] = band

        # ubyte pixel_class(y, x) ;
        #          pixel_class:_FillValue = 0UB ;
        #          pixel_class:long_name = "Pixel classification flags" ;
        #          pixel_class:flag_meanings = "NO_DATA CLEAR_LAND_OR_VEGETATION CLEAR_OCEAN_WATER CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW AMBIGUOUS_CLOUD CLOUD AC_OUT_OF_BOUNDS" ;
        #          pixel_class:flag_values = 0UB, 1UB, 2UB, 3UB, 4UB, 5UB, 6UB, 7UB, 8UB, 9UB ;
        #          pixel_class:coordinates = "crs" ;
        band_name = "pixel_class"
        data = l2w_prototype[band_name] if l2w_prototype else da.random(size=shape, dtype=np.ubyte) % 10
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Pixel classification flags",
                "flag_values": [np.ubyte(0), np.ubyte(1), np.ubyte(2), np.ubyte(3), np.ubyte(4), np.ubyte(5), np.ubyte(6), np.ubyte(7), np.ubyte(8), np.ubyte(9)],
                "flag_meanings": "NO_DATA CLEAR_LAND_OR_VEGETATION CLEAR_OCEAN_WATER CLEAR_INLAND_WATER SNOW_ICE CIRRUS CLOUD_OR_MOUNTAIN_SHADOW AMBIGUOUS_CLOUD CLOUD AC_OUT_OF_BOUNDS",
                "_FillValue": np.ubyte(0),
                "grid_mapping": "crs",
            }
        )
        l1c["quality/water/" + band_name] = band

        # ubyte sen2water_flags(y, x) ;
        #         string sen2water_flags:long_name = "quality and algorithm flags" ;
        #         string sen2water_flags:flag_meanings = "c2rcc_oor acolite_negatives polymer_invalid with_c2rcc with_acolite with_polymer" ;
        #         sen2water_flags:flag_masks = 1UB, 2UB, 4UB, 8UB, 16UB, 32UB ;
        #         string sen2water_flags:coordinates = "crs" ;
        band_name = "sen2water_flags"
        data = l2w_prototype[band_name] if l2w_prototype else da.random(size=shape, dtype=np.ubyte) % 32
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Quality and algorithm flags",
                "flag_masks": [np.ubyte(0), np.ubyte(1), np.ubyte(2), np.ubyte(4), np.ubyte(8), np.ubyte(16), np.ubyte(32)],
                "flag_meanings": "c2rcc_oor acolite_negatives polymer_invalid with_c2rcc with_acolite with_polymer",
                "grid_mapping": "crs",
            }
        )
        l1c["quality/water/" + band_name] = band

        # int pixel_classif_flags(y, x) ;
        #          string pixel_classif_flags:flag_meanings = "IDEPIX_INVALID IDEPIX_CLOUD IDEPIX_CLOUD_AMBIGUOUS IDEPIX_CLOUD_SURE IDEPIX_CLOUD_BUFFER IDEPIX_CLOUD_SHADOW IDEPIX_SNOW_ICE IDEPIX_BRIGHT IDEPIX_WHITE IDEPIX_COASTLINE IDEPIX_LAND IDEPIX_CIRRUS_SURE IDEPIX_CIRRUS_AMBIGUOUS IDEPIX_CLEAR_LAND IDEPIX_CLEAR_WATER IDEPIX_WATER IDEPIX_BRIGHTWHITE IDEPIX_VEG_RISK IDEPIX_MOUNTAIN_SHADOW IDEPIX_POTENTIAL_SHADOW IDEPIX_CLUSTERED_CLOUD_SHADOW" ;
        #          pixel_classif_flags:flag_masks = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 ;
        #          string pixel_classif_flags:flag_coding_name = "pixel_classif_flags" ;
        #          string pixel_classif_flags:flag_descriptions = "Invalid pixels\tPixels which are either cloud_sure or cloud_ambiguous\tSemi transparent clouds, or clouds where the detection level is uncertain\tFully opaque clouds with full confidence of their detection\tA buffer of n pixels around a cloud. n is a user supplied parameter. Applied to pixels masked as \'cloud\'\tPixel is affected by a cloud shadow (combination of shifted cloud mask in cloud gaps and dark clusters coinciding with a corrected shifted cloud mask)\tClear snow/ice pixels\tBright pixels\tWhite pixels\tPixels at a coastline\tLand pixels\tCirrus clouds with full confidence of their detection\tCirrus clouds, or clouds where the detection level is uncertain\tClear land pixels\tClear water pixels\tWater pixels\t\'Brightwhite\' pixels\tPixels with vegetation risk\tPixel is affected by mountain shadow\tPotentially a cloud shadow pixel\tCloud shadow identified by clustering algorithm" ;
        #          pixel_class:coordinates = "crs" ;
        band_name = "pixel_classif_flags"
        data = l2w_prototype[band_name] if l2w_prototype else da.random(size=shape, dtype=np.int)
        band = xr.DataArray(
            data=data,
            dims=("y", "x"),
            attrs={
                "long_name": "Idepix pixel classification flags",
                "flag_values": [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576],
                "flag_meanings": "IDEPIX_INVALID IDEPIX_CLOUD IDEPIX_CLOUD_AMBIGUOUS IDEPIX_CLOUD_SURE IDEPIX_CLOUD_BUFFER IDEPIX_CLOUD_SHADOW IDEPIX_SNOW_ICE IDEPIX_BRIGHT IDEPIX_WHITE IDEPIX_COASTLINE IDEPIX_LAND IDEPIX_CIRRUS_SURE IDEPIX_CIRRUS_AMBIGUOUS IDEPIX_CLEAR_LAND IDEPIX_CLEAR_WATER IDEPIX_WATER IDEPIX_BRIGHTWHITE IDEPIX_VEG_RISK IDEPIX_MOUNTAIN_SHADOW IDEPIX_POTENTIAL_SHADOW IDEPIX_CLUSTERED_CLOUD_SHADOW",
                "_FillValue": np.int(0),
                "grid_mapping": "crs",
            }
        )
        l1c["quality/water/" + band_name] = band

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
        l1c["quality/water/" + band_name] = band

        return inputs
