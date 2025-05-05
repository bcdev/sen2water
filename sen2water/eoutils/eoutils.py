# -*- coding: utf-8 -*-

"""Generic functions to copy Dataset elements"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

from typing import Any

import numpy as np
import xarray as xr


def copy_variable(source: xr.DataArray) -> xr.DataArray:
    """
    Copies a variable into a group of a different product without duplicating data or attributes.

    Parameters
    ----------
    source: EOVariable
        the variable to copy, part of a source EOProduct
    target: EOGroup
        the group to place the new EOVariable into, part of a target EOProduct
    Returns
    -------
    EOVariable
        the variable of the target product just added
    """
    #return source.copy(deep=False, data=source.data)
    return xr.DataArray(source.data, dims=source.dims, attrs=dict(source.attrs))


def copy_attrs(source: Any) -> Any:
    """
    Deep-copies an attribute tree (but not the leaves) and returns the copy.

    Parameters
    ----------
    source: Any
        either a dict or a list or a leaf of the attribute tree
    Returns
    -------
    Any
        a copy of the dict or just a reference to the value.
    """
    return dict(source)


class ScaleAttrs(object):
    def __init__(self, scale_factor, add_offset, fill_value, dtype=np.float32):
        self.scale_factor = scale_factor
        self.add_offset = add_offset
        self.fill_value = fill_value
        self.dtype = dtype

    def scale(self, int_array: np.ndarray):
        float_array = (int_array * self.scale_factor + self.add_offset).astype(
            self.dtype
        )
        if self.fill_value is not None:
            float_array[int_array == self.fill_value] = np.nan
        return float_array

    def unscale(self, float_array: np.ndarray):
        int_array = ((float_array - self.add_offset) / self.scale_factor).astype(
            np.int16
        )
        int_array[np.isnan(float_array)] = (
            self.fill_value if self.fill_value is not None else 0
        )
        return int_array

    @staticmethod
    def create(variable: xr.DataArray):
        if (
            variable.data.dtype != np.float32
            and variable.data.dtype != np.float64
            and ("scale_factor" in variable.attrs or "add_offset" in variable.attrs)
        ):
            return ScaleAttrs(
                (
                    variable.attrs["scale_factor"]
                    if "scale_factor" in variable.attrs
                    else 1.0
                ),
                variable.attrs["add_offset"] if "add_offset" in variable.attrs else 0.0,
                (
                    variable.attrs["_FillValue"]
                    if "_Fill_value" in variable.attrs
                    else None
                ),
            )
        else:
            return None
