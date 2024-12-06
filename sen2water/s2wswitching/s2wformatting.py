# -*- coding: utf-8 -*-

"""Combines Idepix flags and additional flags and selected bands into pixel class"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.51"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 0.51:
# remove statistics attribute, there are individual ones

import uuid
import re
from datetime import datetime
import xarray as xr
from sen2water.eoutils.eoprocessingifc import Operator
from sen2water.eoutils.eoutils import copy_variable
from sen2water.s2wswitching.acolitebands import AcoliteDataset

OCEAN_WAVELENGTHS = [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]
C2RCC_BANDS = [
    "rhow_B1",
    "rhow_B2",
    "rhow_B3",
    "rhow_B4",
    "rhow_B5",
    "rhow_B6",
    "rhow_B7",
    "rhow_B8A",
]
ACOLITE_BANDS = [
    "rhos_443",
    "rhos_492",
    "rhos_560",
    "rhos_665",
    "rhos_704",
    "rhos_740",
    "rhos_783",
    "rhos_833",
    "rhos_865",
    "rhos_1614",
    "rhos_2202",
]
POLYMER_BANDS = [
    "Rw443",
    "Rw490",
    "Rw560",
    "Rw665",
    "Rw705",
    "Rw740",
    "Rw783",
    "Rw842",
    "Rw865",
    "Rw945",
    "Rw1375",
    "Rw1610",
    "Rw2190",
]


class S2wFormatting(Operator):
    """
    Combines output bands, computes statistics and metadata
    """

    def run(
        self,
        water: xr.Dataset,
        idepix: xr.Dataset,
        resampled: xr.Dataset,
        pixelclass: xr.Dataset,
        ocean: xr.Dataset,
        c2rcc: xr.Dataset,
        acolite: AcoliteDataset,
        polymer: xr.Dataset,
        input_id: str,
        with_copyinputs: bool = False,
    ) -> xr.Dataset:
        """
        Creates L2W Dataset with Rw*, pixel_class, sen2water_flags, pixel_classif_flags,
        optionally the individual AC outputs for debugging

        Parameters
        ----------
        water:
            switched reflectances, pixel_class, sen2water_flags
        idepix:
            pixel_classif_flags
        resampled:
            B2, B3, B4, B8A
        pixelclass:
            pixel_class after HR-OC switching
        ocean:
            sen2water_flags after HR-OC switching
        c2rcc:
            C2RCC output with rhow_B1, .., rhow_B8A, c2rcc_flags
        acolite:
            ACOLITE output with rhos_443, .., rhos_2202
        polymer:
            POLYMER output with Rw443, .., Rw2190, bitmask
        input_id:
            S2._MSIL1C_... logical identifier of L1B
        with_copyinputs:
            whether to include individual AC outputs

        Returns
        -------
        xr.Dataset
            Output product with Rw443, .., Rw2190, pixel_class, sen2water_flags,
            pixel_classif_flags,
            optionally rhow_B1_c, .., rhos_443_a, .., Rw443_p, .., c2rcc_flags,
            bitmask, hroc_class, hroc_flags, B2, B3, B4, B8A
        """

        coordinate_bands = {}
        coordinate_bands["y"] = copy_variable(resampled["y"])
        coordinate_bands["x"] = copy_variable(resampled["x"])
        coordinate_bands["crs"] = copy_variable(resampled["crs"])
        # make SNAP happy
        if "wkt" not in coordinate_bands["crs"].attrs:
            coordinate_bands["crs"].attrs["wkt"] = coordinate_bands["crs"].attrs["crs_wkt"]
        # GeoTransform 699960.0 60.0 0.0 5200020.0 0.0 -60.0
        # i2m 60.0,0.0,0.0,-60.0,699960.0,5200020.0
        if "i2m" not in coordinate_bands["crs"].attrs:
            geo_t = coordinate_bands["crs"].attrs["GeoTransform"].split(" ")
            coordinate_bands["crs"].attrs["i2m"] = ",".join([geo_t[i] for i in [1,2,4,5,0,3]])

        target_bands = {}

        for i in range(len(OCEAN_WAVELENGTHS)):
            wavelength = OCEAN_WAVELENGTHS[i]
            bandname = f"Rw{wavelength}"
            target_bands[bandname] = copy_variable(water[bandname])
        target_bands["pixel_class"] = copy_variable(water["pixel_class"])
        target_bands["sen2water_flags"] = copy_variable(water["sen2water_flags"])
        target_bands["pixel_classif_flags"] = copy_variable(
            idepix["pixel_classif_flags"]
        )

        if with_copyinputs:
            target_bands["pixel_class_ac"] = copy_variable(pixelclass["pixel_class"])
            target_bands["sen2water_flags_ac"] = copy_variable(ocean["sen2water_flags"])
            for bandname in C2RCC_BANDS:
                target_bands[f"{bandname}_c"] = copy_variable(c2rcc[bandname])
                target_bands[f"{bandname}_c"].attrs["wavelength_unit"] = "nm"
                del target_bands[f"{bandname}_c"].attrs["valid_pixel_expression"]
            target_bands["c2rcc_flags_c"] = copy_variable(c2rcc["c2rcc_flags"])
            for bandname in acolite.acolite_bands():
                target_bands[f"{bandname}_a"] = copy_variable(acolite[bandname])
                target_bands[f"{bandname}_a"].attrs["radiation_wavelength"] = acolite[
                    bandname
                ].attrs["wavelength"]
                target_bands[f"{bandname}_a"].attrs["wavelength_unit"] = "nm"
            for bandname in POLYMER_BANDS:
                target_bands[f"{bandname}_p"] = copy_variable(polymer[bandname])
                target_bands[f"{bandname}_p"].attrs["radiation_wavelength"] = float(
                    bandname[2:]
                )
                target_bands[f"{bandname}_p"].attrs["wavelength_unit"] = "nm"
            target_bands["polymer_bitmask_p"] = copy_variable(polymer["bitmask"])
            target_bands["B2"] = copy_variable(resampled["B2"])
            target_bands["B3"] = copy_variable(resampled["B3"])
            target_bands["B4"] = copy_variable(resampled["B4"])
            target_bands["B8A"] = copy_variable(resampled["B8A"])

        m = re.match(
            "S2(.)_MSIL1C_(........T......)_N(....)_R(...)_T(.....)_........T......",
            input_id,
        )
        platform = m.group(1)
        start = m.group(2)
        start_date = datetime.strptime(start, "%Y%m%dT%H%M%S")
        start_snapformat = start_date.strftime("%Y-%b-%d %H:%M:%S.%f")
        product_version = m.group(3)
        relorbit = m.group(4)
        granule = m.group(5)
        processing_time = datetime.utcnow().strftime("%Y%m%dT%H%M%S")
        tracking_id = str(uuid.uuid1())
        parameters = "--copyinputs" if with_copyinputs else ""
        processor_version = "0.5"

        target_attrs = {}
        target_attrs["id"] = (
            f"S2{platform}_MSIL2W_{start}_N{product_version}_R{relorbit}_T{granule}_"
            + f"{processing_time}.nc"
        )
        target_attrs["local_id"] = f"T{granule}_{start}_AQU_60m.nc"
        target_attrs["date_created"] = f"{processing_time}Z"
        target_attrs["tracking_id"] = tracking_id
        target_attrs["title"] = "OPT-MPC Sen2Water water-leaving reflectances in 60m"
        target_attrs["institution"] = (
            "Brockmann Consult GmbH as part of OPT-MPC for ESA"
        )
        target_attrs["source"] = f"Sentinel-2 MSI L1C"
        target_attrs["auxiliary"] = "Copernicus 90m DEM,Sen2Water ocean-inland-mask"
        target_attrs["input"] = f"{input_id}"
        target_attrs["processor"] = f"Sen2Water v{processor_version}"
        target_attrs["parameters"] = f"{parameters}"
        target_attrs["product_version"] = f"{product_version[:2]}.{product_version[2:]}"
        target_attrs["history"] = (
            f"Sen2Water switching v{processor_version}; Polymer 4.17beta2; Acolite "
            + f"20241112; C2RCC 9.0cv; Idepix 9.0cv; SNAP-9.0cv; msiresampling v{processor_version}"
        )
        #target_attrs["statistics"] = "TBD"
        target_attrs["references"] = (
            "https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2/"
        )
        target_attrs["license"] = "TBD"
        target_attrs["summary"] = (
            "The Sen2Water L2W product has been processed from Sentinel-2 MSI L1C "
            + "by pixel identification, atmospheric correction with different "
            + "algorithms, and switching between them for ocean and inland water "
            + "pixels. The L2W is provided as part of the Sentinel-2 MSI L2A product "
            + "parallel to the output of Sen2Cor. It can also be generated from L1C "
            + "with Sen2Water stand-alone."
        )
        target_attrs["keywords"] = (
            "reflectance, surface water, ocean optics, Copernicus"
        )
        target_attrs["keywords_vocabulary"] = (
            "NASA Global Change Master Directory (GCMD) Science keywords"
        )
        target_attrs["Conventions"] = "CF-1.11"
        target_attrs["naming_authority"] = "www.esa.int"
        target_attrs["standard_name_vocabulary"] = (
            "NetCDF Climate and Forecast (CF) Metadata Convention"
        )
        target_attrs["creator_name"] = "Brockmann Consult GmbH as part of OPT-MPC for ESA"
        target_attrs["creator_url"] = "https://www.brockmann-consult.de/"
        target_attrs["creator_email"] = "info@brockmann-consult.de"
        target_attrs["contact"] = "TBD"
        target_attrs["project"] = "Optical Mission Performance Centre OPT-MPC"
        target_attrs["cmd_data_type"] = "Grid"
        target_attrs["spatial_resolution"] = "60m"
        target_attrs["time_coverage_start"] = f"{start}Z"
        target_attrs["platform"] = f"Sentinel-2{platform}"
        target_attrs["sensor"] = "MSI"
        target_attrs["start_date"] = start_snapformat
        target_attrs["stop_date"] = start_snapformat
        target_attrs["metadata_profile"] = "beam"
        if with_copyinputs:
            target_attrs["auto_grouping"] = "rhow_B*_c:rhos_*_a:Rw*_p:Rw*:B*"
        return xr.Dataset(target_bands, coordinate_bands, target_attrs)
