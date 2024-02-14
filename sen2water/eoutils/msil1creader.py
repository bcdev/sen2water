# -*- coding: utf-8 -*-

"""MSI L1C reader implementation as entrypoint (engine) for xarray"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import re
import glob
import os.path
import warnings
from datetime import datetime
import operator
from functools import reduce
from lxml import etree
import xarray as xr
import numpy as np
from xarray.backends import BackendEntrypoint


class MsiL1cBackendEntrypoint(BackendEntrypoint):

    band_name_in_file = {
        "B01": "B1",
        "B02": "B2",
        "B03": "B3",
        "B04": "B4",
        "B05": "B5",
        "B06": "B6",
        "B07": "B7",
        "B08": "B8",
        "B8A": "B8A",
        "B09": "B9",
        "B10": "B10",
        "B11": "B11",
        "B12": "B12",
    }
    band_names = [
        "B1",
        "B2",
        "B3",
        "B4",
        "B5",
        "B6",
        "B7",
        "B8",
        "B8A",
        "B9",
        "B10",
        "B11",
        "B12",
    ]
    resolutions = {
        "B1": 60,
        "B2": 10,
        "B3": 10,
        "B4": 10,
        "B5": 20,
        "B6": 20,
        "B7": 20,
        "B8": 10,
        "B8A": 20,
        "B9": 60,
        "B10": 60,
        "B11": 20,
        "B12": 20,
    }
    wavelengths = {
        "B1": 443.0,
        "B2": 490.0,
        "B3": 560.0,
        "B4": 665.0,
        "B5": 705.0,
        "B6": 740.0,
        "B7": 783.0,
        "B8": 842.0,
        "B8A": 865.0,
        "B9": 945.0,
        "B10": 1375.0,
        "B11": 1610.0,
        "B12": 2190.0,
    }
    quality_flags = [
        "ancillary_lost",
        "ancillary_degraded",
        "msi_lost",
        "msi_degraded",
        "defective",
        "nodata",
        "partially_corrected_crosstalk",
        "saturated_l1a",
    ]
    cloud_flags = [
        "opaque_clouds",
        "cirrus_clouds",
        "snow_and_ice_areas"
    ]

    var_pattern = ".*_(...).jp2"

    chunksize_in_meters = 36600

    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        merge_flags=False,
    ):
        # TODO make chunking a parameter and consider it here and in resampling during writing
        # TODO move variable attributes from resampling to here, reuse them in resampling instead
        # TODO parse metadata file and convert into global attributes

        datasets = []
        granule_dir = glob.glob(os.path.join(filename_or_obj, "GRANULE/*"))[0]

        # read metadata from XML metadata
        mtd_path = os.path.join(filename_or_obj, "MTD_MSIL1C.xml")
        l1c_dom = etree.parse(mtd_path)
        scale_factor = 1.0 / int(
            l1c_dom.xpath(
                "n1:General_Info/Product_Image_Characteristics/QUANTIFICATION_VALUE",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
                },
            )[0].text
        )
        add_offsets = [int(offset.text) * scale_factor for offset in l1c_dom.xpath(
            "n1:General_Info/Product_Image_Characteristics/Radiometric_Offset_List/RADIO_ADD_OFFSET",
            namespaces={
                "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
            }
        )]
        #wavelengths = [np.round(float(w.text)) for w in l1c_dom.xpath(
        #    "n1:General_Info/Product_Image_Characteristics/Spectral_Information_List/Spectral_Information/Wavelength/CENTRAL",
        #    namespaces={
        #        "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
        #    }
        #)]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            # open reflectances
            for refl_path in glob.glob(os.path.join(granule_dir, "IMG_DATA/*B??.jp2")):
                filename = refl_path[refl_path.rfind("/") + 1 :]
                variable = self.band_name_in_file[re.match(self.var_pattern, filename).group(1)]
                band_index = self.band_names.index(variable)
                resolution = self.resolutions[variable]
                chunks = self.chunksize_in_meters // resolution
                #print(f"{variable=} {chunks=}")
                self.open_reflectance_bands(refl_path, band_index, variable, resolution, chunks, scale_factor, add_offsets,
                                            datasets)
            # open detector footprint
            for refl_path in glob.glob(
                os.path.join(granule_dir, "QI_DATA/MSK_DETFOO_*.jp2")
            ):
                filename = refl_path[refl_path.rfind("/") + 1 :]
                variable = self.band_name_in_file[re.match(self.var_pattern, filename).group(1)]
                resolution = self.resolutions[variable]
                chunks = self.chunksize_in_meters // resolution
                self.open_detector_footprints(refl_path, variable, resolution, chunks, datasets)
            # open quality masks
            for refl_path in glob.glob(
                os.path.join(granule_dir, "QI_DATA/MSK_QUALIT_*.jp2")
            ):
                filename = refl_path[refl_path.rfind("/") + 1 :]
                variable = self.band_name_in_file[re.match(self.var_pattern, filename).group(1)]
                resolution = self.resolutions[variable]
                chunks = self.chunksize_in_meters // resolution
                self.open_quality_flags(refl_path, variable, resolution, chunks, merge_flags, datasets)
            # open cloud masks
            for refl_path in glob.glob(
                os.path.join(granule_dir, "QI_DATA/MSK_CLASSI_B00.jp2")
            ):
                resolution = self.resolutions["B1"]
                chunks = self.chunksize_in_meters // resolution
                self.open_cloud_masks(refl_path, resolution, chunks, merge_flags, datasets)
            # open atmospheric ancillary data
            for auxfile in ["AUX_CAMSFO", "AUX_ECMWFT"]:
                auxpath = os.path.join(granule_dir, f"AUX_DATA/{auxfile}")
                self.open_ancillary_files(auxpath, datasets)
            # read angles from XML metadata
            mtd_path = os.path.join(granule_dir, "MTD_TL.xml")
            self.open_angles(mtd_path, datasets)

            vars = reduce(operator.ior, [d.data_vars for d in datasets], {})
            attrs = reduce(operator.ior, [d.attrs for d in datasets], {})
            start_time = datetime.strptime(l1c_dom.xpath(
                "n1:General_Info/Product_Info/PRODUCT_START_TIME",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
                },
            )[0].text[:-1]+"000", "%Y-%m-%dT%H:%M:%S.%f")
            stop_time = datetime.strptime(l1c_dom.xpath(
                "n1:General_Info/Product_Info/PRODUCT_STOP_TIME",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
                },
            )[0].text[:-1]+"000", "%Y-%m-%dT%H:%M:%S.%f")
            spacecraft = l1c_dom.xpath(
                "n1:General_Info/Product_Info/Datatake/SPACECRAFT_NAME",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd"
                },
            )[0].text
            attrs["start_date"] = start_time.strftime("%d-%b-%Y %H:%M:%S.%f").upper()
            attrs["stop_date"] = stop_time.strftime("%d-%b-%Y %H:%M:%S.%f").upper()
            attrs["spacecraft"] = spacecraft

            ds = xr.Dataset(data_vars=vars, attrs=attrs)
            return ds

    def open_angles(self, mtd_path, datasets):
        dom = etree.parse(mtd_path)
        sza = self.get_values(
            dom,
            "n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid/Zenith/Values_List/VALUES",
        )
        saa = self.get_values(
            dom,
            "n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid/Azimuth/Values_List/VALUES",
        )
        bands = {
            "sun_zenith": xr.DataArray(sza, dims=["y_tp", "x_tp"]),
            "sun_azimuth": xr.DataArray(saa, dims=["y_tp", "x_tp"])
        }
        for band_id in range(13):
            vza_data = []  # np.empty((12, 23, 23))
            vaa_data = []  # np.empty((12, 23, 23))
            # vza_data[:] = np.nan
            # vaa_data[:] = np.nan
            detectors = []
            for detector in range(1, 13):
                detector_vza = self.get_values(
                    dom,
                    'n1:Geometric_Info/Tile_Angles/Viewing_Incidence_Angles_Grids[@bandId="{}" and @detectorId="{}"]/Zenith/Values_List/VALUES'.format(
                        band_id, detector
                    ),
                )
                detector_vaa = self.get_values(
                    dom,
                    'n1:Geometric_Info/Tile_Angles/Viewing_Incidence_Angles_Grids[@bandId="{}" and @detectorId="{}"]/Azimuth/Values_List/VALUES'.format(
                        band_id, detector
                    ),
                )
                if detector_vza is not None and detector_vaa is not None:
                    vza_data.append(detector_vza)
                    vaa_data.append(detector_vaa)
                    detectors.append(detector)
            bands["vza_{}".format(self.band_names[band_id])] = xr.DataArray(
                np.stack(vza_data), dims=["detector", "y_tp", "x_tp"]
            )
            bands["vaa_{}".format(self.band_names[band_id])] = xr.DataArray(
                np.stack(vaa_data), dims=["detector", "y_tp", "x_tp"]
            )
            bands["detector"] = xr.DataArray(
                np.array(detectors, dtype=np.uint8), dims=["detector"]
            )
        datasets.append(xr.Dataset(bands))
        shape_y_x = self.get_shape(
            dom,
            "n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid/Zenith/Values_List/VALUES",
        )
        ymax = float(
            dom.xpath(
                'n1:Geometric_Info/Tile_Geocoding/Geoposition[@resolution="10"]/ULY',
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
                },
            )[0].text
        )
        xmin = float(
            dom.xpath(
                'n1:Geometric_Info/Tile_Geocoding/Geoposition[@resolution="10"]/ULX',
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
                },
            )[0].text
        )
        ystep = float(
            dom.xpath(
                "n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid/Zenith/ROW_STEP",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
                },
            )[0].text
        )
        xstep = float(
            dom.xpath(
                "n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid/Zenith/COL_STEP",
                namespaces={
                    "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
                },
            )[0].text
        )
        # y = [ymax - i * ystep - ystep / 2 for i in range(shape_y_x[0])]
        # x = [xmin + i * xstep + xstep / 2 for i in range(shape_y_x[1])]
        y = [ymax - i * ystep for i in range(shape_y_x[0])]
        x = [xmin + i * xstep for i in range(shape_y_x[1])]
        datasets.append(xr.Dataset({"y_tp": y, "x_tp": x}))

    def open_ancillary_files(self, auxpath, datasets):
        ds = xr.open_dataset(auxpath, engine="cfgrib", mask_and_scale=False)
        ds = ds.rename_dims(
            {"latitude": "aux_latitude", "longitude": "aux_longitude"}
        ).rename_vars({"latitude": "aux_latitude", "longitude": "aux_longitude"})
        datasets.append(ds)

    def open_cloud_masks(self, refl_path, resolution, chunks, merge_flags, datasets):
        ds = xr.open_dataset(
            refl_path,
            chunks={"y": chunks, "x": chunks},
            engine="rasterio",
            mask_and_scale=False,
        )
        ds = ds.rename_dims(
            {"y": f"y{resolution}", "x": f"x{resolution}"}
        ).rename_vars({"y": f"y{resolution}", "x": f"x{resolution}"})
        da = ds["band_data"].drop_vars("band")
        if merge_flags:
            target_variable = "cloud_ice_flags"
            attrs = da.attrs
            attrs["flag_masks"] = [np.int32(1), np.int32(2), np.int32(4)]
            attrs["flag_meanings"] = [flag for flag in self.cloud_flags]
            datasets.append(
                xr.Dataset(
                    {
                        target_variable: xr.DataArray(
                            da.data[0] + 2 * da.data[1] + 4 * da.data[2],
                            da.coords,
                            da.dims[1:],
                            target_variable,
                            attrs,
                        )
                    }
                )
            )
        else:
            for i, flag in enumerate(self.cloud_flags):
                datasets.append(
                    xr.Dataset(
                        {
                            f"B_{flag}": xr.DataArray(
                                da.data[i],
                                da.coords,
                                da.dims[1:],
                                flag,
                                da.attrs,
                            )
                        }
                    )
                )
        del da, ds

    def open_reflectance_bands(self, refl_path, band_index, variable, resolution, chunks, scale_factor, add_offsets,
                               datasets):
        ds = xr.open_dataset(
            refl_path,
            chunks={"y": chunks, "x": chunks},
            engine="rasterio",
            mask_and_scale=False,
        )
        ds = ds.rename_dims(
            {"y": f"y{resolution}", "x": f"x{resolution}"}
        ).rename_vars({"y": f"y{resolution}", "x": f"x{resolution}"})
        da = ds["band_data"].drop_vars("band")
        attrs = {
            "long_name": f"Reflectance in band {variable}",
            "units": "dl",
            "scale_factor": scale_factor,
            "add_offset": add_offsets[band_index],
            "radiation_wavelength": self.wavelengths[variable],
            "radiation_wavelength_unit": "nm",
            "spectral_band_index": np.int32(band_index)
        }
        datasets.append(
            xr.Dataset(
                {
                    variable: xr.DataArray(
                        da.data[0], da.coords, da.dims[1:], variable, attrs
                    )
                }
            )
        )
        del da, ds

    def open_detector_footprints(self, refl_path, variable, resolution, chunks, datasets):
        ds = xr.open_dataset(
            refl_path,
            chunks={"y": chunks, "x": chunks},
            engine="rasterio",
            mask_and_scale=False,
        )
        ds = ds.rename_dims(
            {"y": f"y{resolution}", "x": f"x{resolution}"}
        ).rename_vars({"y": f"y{resolution}", "x": f"x{resolution}"})
        da = ds["band_data"].drop_vars("band")
        target_variable = f"B_detector_footprint_{variable}"
        datasets.append(
            xr.Dataset(
                {
                    target_variable: xr.DataArray(
                        da.data[0],
                        da.coords,
                        da.dims[1:],
                        target_variable,
                        da.attrs,
                    )
                }
            )
        )
        del da, ds

    def open_quality_flags(self, refl_path, variable, resolution, chunks, merge_flags, datasets):
        ds = xr.open_dataset(
            refl_path,
            chunks={"y": chunks, "x": chunks},
            engine="rasterio",
            mask_and_scale=False,
        )
        ds = ds.rename_dims(
            {"y": f"y{resolution}", "x": f"x{resolution}"}
        ).rename_vars({"y": f"y{resolution}", "x": f"x{resolution}"})
        da = ds["band_data"].drop_vars("band")
        if merge_flags:
            target_variable = f"quality_flags_{variable}"
            attrs = da.attrs
            attrs["flag_masks"] = [np.int32(1), np.int32(2), np.int32(4), np.int32(8), np.int32(16), np.int32(32), np.int32(64), np.int32(128)]
            attrs["flag_meanings"] = [flag for flag in self.quality_flags]
            datasets.append(
                xr.Dataset(
                    {
                        target_variable: xr.DataArray(
                            (da.data[0] + 2 * da.data[1] + 4 * da.data[2] + 8 * da.data[3] + 16 * da.data[4] +
                             32 * da.data[5] + 64 * da.data[6] + 128 * da.data[7]).astype(np.uint8),
                            da.coords,
                            da.dims[1:],
                            target_variable,
                            attrs,
                        )
                    }
                )
            )
        else:
            for i, prefix in enumerate(self.quality_flags):
                target_variable = f"B_{prefix}_{variable}"
                datasets.append(
                    xr.Dataset(
                        {
                            target_variable: xr.DataArray(
                                da.data[i],
                                da.coords,
                                da.dims[1:],
                                target_variable,
                                da.attrs,
                            )
                        }
                    )
                )
        del da, ds

    open_dataset_parameters = ["filename_or_obj", "drop_variables"]

    def guess_can_open(self, filename_or_obj):
        try:
            return os.exists(os.path.join(filename_or_obj, "MTD_MSIL1C.xml"))
        except TypeError:
            return False

    description = "Read MSI L1C data products in SAFE format"

    url = "https://link_to/your_backend/documentation"

    @staticmethod
    def get_values(dom, xpath):
        list = dom.xpath(
            xpath,
            namespaces={
                "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
            },
        )
        if len(list) > 0:
            return np.array([[np.float32(i) for i in x.text.split()] for x in list], dtype=np.float32)
        else:
            return None

    @staticmethod
    def get_shape(dom, xpath):
        list = dom.xpath(
            xpath,
            namespaces={
                "n1": "https://psd-14.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd"
            },
        )
        return [len(list), len(list[0].text.split())]
