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

from typing import Any, Optional, Literal, Tuple, Dict, Iterable, List
import numpy as np
import dask.array as da
import xarray as xr

from pyproj import Transformer, CRS
from eopf.computing import EOProcessingUnit, MappingDataType, MappingAuxiliary

from sen2water.msiidepix.libs.algorithms import CloudMasks
from sen2water.msiresampling.constants import MsiConstantsReengineering, MsiConstants


class IdepixClassificationPU(EOProcessingUnit):

    def run(
            self,
            inputs: MappingDataType,
            adfs: Optional[MappingAuxiliary] = None,
            resolution: int = 60,
            chunksize_in_meters=36600,
            **kwargs,
    ) -> MappingDataType:

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c: xr.DataTree = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        input_data_list = []


        # determine parameters
        dims = {
            "y": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["x"]
        }
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)

        # create result accumulator
        pixel_classif_flag = xr.DataTree(name="pixel_classif_flag")

        # add geo-coding with 1-d metric coordinates, CRS, 2-d geographic coordinates

        # geographic_coords_dt = self._add_geocoding(
        #     l1c,
        #     resampled,
        #     ancillary,
        #     resolution,
        #     dims,
        #     chunks,
        # )

        # retrieval
        pixel_classif = da.map_blocks(
            self.compute_pixel_classification_flag,
            l1c,
            dtype=dtype,
            meta=np.array((), dtype=dtype),
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
            l1c: xr.DataTree,
            *,
            image_shape: Tuple[int, int],
            image_chunksize: Tuple[int, int],
            block_id: Tuple[int, int],
    ) -> np.ndarray:
        """
        Block-wise interpolation of tie-point data to the image grid
        using linear interpolation based on pixel coordinates.

        Parameters
        ----------
        l1c: xr.DataTree()
            full source data

        Returns
        -------
        np.ndarray
            array of pixel classif flags with the extent of the image block,
            which is usually image_chunksize except for the right and lower border.
        """

        refl_bands = MsiConstants.bands

        for i in range(len(refl_bands)):
            refl = l1c[refl_bands[i]].data

        dum = da.full_like(qua, False)  # dummy
        cloud_masks = CloudMasks(
            dtype=self.dtype, neural_network_path=self.neural_network_path
        ).apply_to(inv, lnd, gli, dum, dum, dum, *toa)



        return result

    func = compute_pixel_classification_flag





