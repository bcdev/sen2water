# -*- coding: utf-8 -*-

"""Chain of operators that switch between C2RCC, ACOLITE, POLYMER and format the L2W"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

from datetime import datetime
from typing import Any, Dict
import xarray as xr
from sen2water.eoutils.eologging import logger
from sen2water.eoutils.eoprocessingifc import Operator
from sen2water.s2wswitching.pixelclassification import PixelClassification
from sen2water.s2wswitching.coastalswitching import CoastalWaterSwitching
from sen2water.s2wswitching.inlandwaterswitching import InlandWaterSwitching
from sen2water.s2wswitching.s2wformatting import S2wFormatting


class SwitchingProcessor(Operator):
    """
    Implementation of Sen2Water AC switching and output formatting processor.
    """

    def run(
        self,
        resampled: xr.Dataset,
        idepix: xr.Dataset,
        c2rcc: xr.Dataset,
        acolite: xr.Dataset,
        polymer: xr.Dataset,
        s2wmask: xr.Dataset,
        input_id: str,
        with_copyinputs: bool = False,
    ) -> xr.Dataset:
        """
        Chains pixel classification, coastal switching, inland water switching, and
        formatting.

        Parameters
        ----------
        resampled:
            60m MSI L1C with B1, .., B12, lat, lon
        idepix:
            Idepix output with pixel_classif_flags
        c2rcc:
            C2RCC output with rhow_B1, .., rhow_B8A, c2rcc_flags
        acolite: str
            ACOLITE output with rhos_443, .., rhos_2202
        polymer: str
            POLYMER output with Rw443, .., Rw2180, bitmask
        s2wmask: str
            Static Sen2Water ocean/inland water/land mask
        with_copyinputs: bool
            whether to add individual AC results and coastal switching flags to the output

        Returns
        -------
        xr.Dataset:
            output product

        """

        # Pixel classification determines pixel_class mainly from pixel_classif_flags and s2wmask
        pixelclass = PixelClassification().run(
            resampled, idepix, c2rcc, acolite, s2wmask
        )
        logger.info("pixel classification prepared")
        # Coastal switching determines sen2water_flags and ocean reflectances switched between c2rcc and acolite
        ocean = CoastalWaterSwitching().run(c2rcc, acolite, pixelclass)
        logger.info("ocean and coastal water switching prepared")
        # InlandWaterSwitching updates sen2water_flags and pixel_class and switches between coastal and polymer
        water = InlandWaterSwitching().run(ocean, polymer, pixelclass, s2wmask)
        logger.info("inland water switching prepared")
        # S2wFormatting assembles reflectances, pixel_class, sen2water_flags, pixel_classif_flags, optionally more
        output = S2wFormatting().run(
            water,
            idepix,
            resampled,
            pixelclass,
            ocean,
            c2rcc,
            acolite,
            polymer,
            input_id,
            with_copyinputs,
        )
        logger.info("output formatting prepared")
        return output
