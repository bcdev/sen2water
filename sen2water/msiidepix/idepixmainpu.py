from typing import Any, Optional, List

import xarray as xr
from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from sen2water.msiidepix.constants import IdepixMsiConstants as ic
from sen2water.msiidepix.idepixclassificationpu import IdepixClassificationPU


class IdepixMainPU(EOProcessingUnit):
    """
    TODO add doc
    """
    PROCESSOR_NAME = "MsiIdepix"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    bands = ic.bands
    resolutions = ic.resolutions

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

        idepix = xr.DataTree(name="idepix")

        # compute Idepix classification...
        self._compute_idepix_classification(l1c, idepix)

        # compute slope, aspect, orientation...
        self._compute_slope_aspect_orientatiomn(l1c, idepix)

        return {
            "idepix": idepix
        }


    def _compute_idepix_classification(
            self,
            l1c: xr.DataTree,
            idepix: xr.DataTree,
    ):
        """
        Computes Idepix 'raw' classification flag from spectral tests
        """

        flag_band_name = "quality/mask/pixel_classif_flag"

        # TEST: simplified Idepix version, works with 'PU --> PU' approach
        # idepix_flags = IdepixClassificationPUPU().run(
        #     {'l1c': l1c},
        # )["pixel_classif_flag"]

        # 'PU --> Algorithm' approach, this is what we finally want
        idepix_flags = IdepixClassificationPU().run(
            {'l1c': l1c},
        )["pixel_classif_flag"]

        idepix[flag_band_name] = self._with_encoding(xr.DataArray(
            idepix_flags["pixel_classif_flag"]
        ))


    def _compute_slope_aspect_orientatiomn(self, l1c, idepix):
        """
        Computes slope aspect orientatiomn from spectral tests
        TODO

        Parameters
        ----------
        l1c
        idepix

        Returns
        -------

        """
        pass


    def _with_encoding(self, var: xr.DataArray) -> xr.DataArray:
        var.encoding = {
            "zlib": True,
            "complevel": 5,
            "chunksizes": var.data.chunksize,
        }
        return var
