from typing import Any, Optional, List

import kwargs
import numpy as np

import xarray as xr
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from msic2rcc.c2rccmsialgorithm import C2rccMsiAlgorithm
from msic2rcc.constants import C2rccMsiConstants as cc


class C2rccMsiPU(EOProcessingUnit):
    """
    TODO add doc
    """
    PROCESSOR_NAME = "MsiC2rcc"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    bands = cc.bands
    resolutions = cc.resolutions

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
            **kwargs
    ) -> MappingDataType:
        """
        TODO: fill doc

        Parameters
        ----------
        inputs
        adfs
        mode
        resolution
        chunksize_in_meters
        ancillary
        kwargs

        Returns
        -------

        """

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c: xr.DataTree = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        # determine parameters if needed...

        # create result accumulator

        c2rcc = xr.DataTree(name="c2rcc")

        # compute C2RCC...
        self._compute_c2rcc(l1c, c2rcc)

        return {
            "c2rcc": c2rcc
        }


    def _compute_c2rcc(
            self,
            l1c: xr.DataTree,
            c2rcc: xr.DataTree,
    ):
        """
        Computes C2RCC quantities
        """
        dims = {
            "y": l1c[f"measurements/reflectance/resampled"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/resampled"].coords.sizes["x"]
        }

        toa_data = [l1c[f"measurements/reflectance/resampled/{band}"].data for band in cc.bands]
        c2rcc_result = C2rccMsiAlgorithm().apply(
            *toa_data,
            thresh_cw=0.007,
            thresh_gcl=-0.11,
            thresh_cl=0.007,
            dtype=np.int32
        )

        c2rcc_quantities = xr.DataTree()
        # TODO: continue


    def _with_encoding(self, var: xr.DataArray) -> xr.DataArray:
        var.encoding = {
            "zlib": True,
            "complevel": 5,
            "chunksizes": var.data.chunksize,
        }
        return var
