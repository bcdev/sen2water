# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Olaf Danne (Brockmann Consult GmbH)"
__copyright__ = "Copyright 2026, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

from typing import Optional, Tuple, Dict, Iterable

import dask.array as da
import numpy as np
import xarray as xr
from eopf.computing import EOProcessingUnit, MappingDataType, MappingAuxiliary

from sen2water.msiidepix.constants import IdepixMsiConstants as ic
from sen2water.msiidepix.libs.msialgorithms import MsiPixelClassification


class IdepixClassificationAlgoPU(EOProcessingUnit):
    def run(
            self,
            inputs: MappingDataType,
            adfs: Optional[MappingAuxiliary] = None,
            dims: Tuple[int, int] = None,
            dtype: str = None,
            image_chunksize: Tuple[int, int] = None,
            **kwargs,
    ) -> MappingDataType:

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        toa_data = [l1c[f"measurements/reflectance/resampled/{band}"].data for band in ic.bands]

        # retrieval
        pixel_classif = da.map_blocks(
            self.compute_pixel_classification_flag,
            *toa_data,
            dtype=np.int32,
            meta=np.array((), dtype=np.int32)
            **kwargs,
        )
        result = xr.DataTree()
        result["pixel_classif_flag"] = xr.DataArray(
            pixel_classif,
            dims=dims,
            attrs={"long_name": "pixel_classif_flag"},
        )
        return {"pixel_classif_flag": result}


    def compute_pixel_classification_flag(
            self,
            *toa: xr.DataTree,
    ) -> np.ndarray:
        """
        Computation of Idepix classification flag

        Parameters
        ----------
        l1c: xr.DataTree()
            source data

        image_chunksize: Tuple[int, int]
            chunk size

        Returns
        -------
        result - np.ndarray
            array of pixel classif flags with the extent of the image block, which is usually image_chunksize
        """

        # classification logic
        tc4_cirrus_value = -0.8239 * toa[1] + 0.0849 * toa[2] + 0.4396 * toa[3] - 0.058 * toa[8] + 0.2013 * toa[
            11] - 0.2773 * toa[12] - toa[10]

        tc4_value = -0.8239 * toa[1] + 0.0849 * toa[2] + 0.4396 * toa[3] - 0.058 * toa[8] + 0.2013 * toa[
            11] - 0.2773 * toa[12]

        ndwi_value = (toa[8] - toa[11]) / (toa[8] + toa[11])

        b3_b11_value = (toa[2] / toa[11])

        vis_bright_value = (toa[1] + toa[2] + toa[3]) / 3.0

        is_b3_b11 = b3_b11_value > ic.B3B11_THRESH

        is_invalid = np.zeros( toa[0].shape, dtype=bool)  # TODO implement
        is_clear_snow = np.zeros(toa[0].shape, dtype=bool)  # TODO implement

        gcw = tc4_cirrus_value < ic.GCW_THRESH
        tcw = self._and(tc4_value < ic.TCW_TC_THRESH, ndwi_value < ic.TCW_NDWI_THRESH)
        acw = self._and(is_b3_b11, (self._or(gcw, tcw)))
        gcl = self._and(self._not(is_b3_b11), tc4_cirrus_value < ic.GCL_THRESH_DEFAULT,
                        vis_bright_value > ic.VISBRIGHT_THRESH)
        is_cloud_sure = self._and(self._not(is_invalid), self._not(is_clear_snow), self._or(acw, gcl))

        is_cloud_ambiguous = b3_b11_value > ic.B3B11_THRESH

        is_cloud = self._or(is_cloud_sure, is_cloud_ambiguous)

        # final classification
        result = np.zeros(shape=toa[0].shape, dtype=np.int32)

        result[is_cloud] = ic.FLAG_MASK_CLOUD
        result[is_cloud_sure] = ic.FLAG_MASK_CLOUD | ic.FLAG_MASK_CLOUD_SURE
        result[is_cloud_ambiguous] = ic.FLAG_MASK_CLOUD | ic.FLAG_MASK_CLOUD_AMBIGUOUS

        return result


    # def _compute_pixel_classification_flag(
    #         self,
    #         *toa: xr.DataTree,
    #         result: xr.DataTree,
    #         ancillary: Iterable[str],
    #         lat_data: da.Array,
    #         lon_data: da.Array,
    #         dims: Dict[str, int],
    #         chunks: Tuple[int,int],
    # ):
    #     """
    #     'Algorithm' approach
    #
    #     :param l1c:
    #     :param resampled:
    #     :param ancillary:
    #     :param lat_data:
    #     :param lon_data:
    #     :param dims:
    #     :param chunks:
    #     :return:
    #     """
    #     pixel_classif_band = l1c[anc_band_name]
    #     pixel_classif_flag_data = MsiPixelClassification().apply(
    #         *toa,
    #         lat_data,
    #         lon_data,
    #     )
    #
    #     anc_lat = l1c["conditions/meteorology/ecmwf/latitude"]
    #     anc_lon = l1c["conditions/meteorology/ecmwf/longitude"]
    #     anc_band = l1c[anc_band_name]
    #     anc_data = AncillaryInterpolation().apply(
    #         lat_data,
    #         lon_data,
    #         anc_lat=anc_lat.values,
    #         anc_lon=anc_lon.values,
    #         anc_data=anc_band.values,
    #         variable=anc_band,
    #         dtype=np.float32,
    #     )
    #     resampled[f"conditions/meteorology/resampled/{anc_band.name}"] = xr.DataArray(
    #         anc_data,
    #         dims=dims,
    #         attrs=anc_band.attrs.copy() | {
    #             "coordinates": "y x",
    #             "grid_mapping": "crs",
    #         }
    #     )
    #
    # func = _compute_pixel_classification_flag