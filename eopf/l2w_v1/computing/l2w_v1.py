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
            raise MyError("Missing mandatory input product for the processor run.")


        # introduce group
        water_leaving_reflectance = xr.DataTree(name="water_leaving_reflectance")
        # print(inputs)
        l1c = inputs["l1c"]
        l1c["measurements/water_leaving_reflectance"] = water_leaving_reflectance
        shape = l1c["measurements/reflectance/r60m/b01"].data.shape

        # introduce variables
        for band_name in ["Rw443", "Rw560", "Rw665"]:
            data = da.random.random(size=shape)
            band = xr.DataArray(data=data, dims=("y", "x"))
            # water_leaving_reflectance[band_name] = band
            l1c["measurements/water_leaving_reflectance/" + band_name] = band

        return inputs
