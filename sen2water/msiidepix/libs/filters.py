#  Copyright (c) 2023 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2023 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

import numpy as np

# noinspection PyPackageRequirements
import scipy as scp

# noinspection PyPackageRequirements
from numba import carray

# noinspection PyPackageRequirements
from numba import cfunc

# noinspection PyPackageRequirements
from numba.types import CPointer

# noinspection PyPackageRequirements
from numba.types import float64

# noinspection PyPackageRequirements
from numba.types import intc

# noinspection PyPackageRequirements
from numba.types import intp

# noinspection PyPackageRequirements
from numba.types import voidptr


@cfunc(intc(CPointer(float64), intp, CPointer(float64), voidptr))
def _nb_max(values_ptr, len_values, result, data):
    """Computes the maximum valid value."""
    values = carray(values_ptr, (len_values,), dtype=float64)
    result[0] = np.nanmax(values)
    return 1


filter_max = scp.LowLevelCallable(_nb_max.ctypes)
"""The maximum filter function."""
