# -*- coding: utf-8 -*-
#  Copyright (c) 2024 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2024 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

import os
import re
import time
from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from typing import Iterable
from typing import Mapping

import netCDF4
import numpy as np

from idepix_msi.eoutils.eoutils import copy_variable


@dataclass
class VariableData:
    values: np.ndarray
    attrs: Mapping
    dims: Iterable


class DatasetAdapter(ABC):
    """
    A wrapper for the ``Dataset`` types of xarray and netCDF4.
    The ``DatasetAdapter`` provides a uniform interface to access variables and their data.
    The ``DatasetAdapter`` class is never instantiated directly. Instead, an instance of
    a subclass for the corresponding dataset is created.
    """

    # TODO this contains only the variables required by IdePix
    group_dict = {
        "altitude": "geo_coordinates.nc",
        "latitude": "geo_coordinates.nc",
        "longitude": "geo_coordinates.nc",
        "detector_index": "instrument_data.nc",
        "solar_flux": "instrument_data.nc",
        "quality_flags": "qualityFlags.nc",
        "SZA": "tie_geometries.nc",
        "SAA": "tie_geometries.nc",
        "OZA": "tie_geometries.nc",
        "OAA": "tie_geometries.nc",
    }

    def __new__(cls, *args, **kwargs):
        ds = args[0]
        if is_netcdf4_type(ds):
            import netCDF4

            return super().__new__(NetCDF4DatasetAdapter)
        elif is_xarray_type(ds):
            import xarray as xr

            return super().__new__(XarrayDatasetAdapter)
        else:
            raise ValueError(f"Invalid dataset type: {repr(type(ds))}")

    @abstractmethod
    def get_variable_data(self, item):
        """
        Data access to retrieve the data of a variable from the
        dataset. Must be implemented by each subclass.
        The returned data must implement the numpy array interface
        and may be e.g. a dask array
        """
        ...

    @abstractmethod
    def get_variable_values(self, item) -> np.ndarray:
        """
        Data access to retrieve the values (np.ndarray) of a variable from the
        dataset. Must be implemented by each subclass.
        The values must implement the numpy array interface and represent concrete values.
        """
        ...

    @abstractmethod
    def get_variable(self, item):
        """
        Returns a ``variable`` type object of the wrapped ``Dataset`` type.
        """
        ...

    @classmethod
    @abstractmethod
    def create_dataset_from_parts(
        cls,
        *,
        variables: Mapping[str, VariableData],
        coordinates: Mapping,
        attributes: Mapping,
        base_ds=None,
    ):
        """
        Generate a new dataset of the wrapped dataset type from
        variables, coordinates and attributes.
        This method provides a uniform interface to the constructors of dataset
        types.
        Must be implemented by each subclass.

        :param variables: Dictionary of variable names mapped
            to ``VariableData`` objects specifying the data for the desired variables
        :type variables: Mapping

        :param coordinates: Dictionary of coordinate names mapped to coordinate
            objects or variables specifying the new coordinates
        :type coordinates: Mapping

        :param attributes: Mapping of attribute names to attribute values.
        :type attributes: Mapping

        :param base_ds: Some types of Dataset (e.g. netCDF4.Dataset) require a ``Dataset`` instance
            to be created before variables are added. If this ``base_ds`` is not ``None``,
            variables and coordinates will be added to ``base_ds`` instead of a new object.
            This is not supported by all ``Dataset`` types.
        :type base_ds: optional
        """
        ...

    @abstractmethod
    def copy_variable(self, variable: str):
        """
        Create a copy of the variable with name ``variable``.
        """
        ...


class NetCDF4DatasetAdapter(DatasetAdapter):
    radiance_re = re.compile("Oa(0[1-9]|1[0-9]|2[0-1])_radiance")

    def __init__(self, ds):
        self.ds = ds

    def get_variable(self, item, group=None) -> netCDF4.Variable:
        group = group or self._find_group(item)
        if group in self.ds.groups.keys():
            return self.ds.groups[group][item]

        if (variable := self.ds.variables.get(item, None)) is not None:
            return variable

        raise KeyError(
            f"Could not find variable {item} in dataset."
            " Currently only the variables required by IdePix are searched."
        )

    def _find_group(self, var: str) -> str | None:
        """
        Returns the group in the netCDF4 dataset that contains the variable
        """
        if self.radiance_re.match(var) is not None:
            return var + ".nc"

        return self.group_dict.get(var)

    def get_variable_data(self, item, *, allow_masked=False):
        var = self.get_variable(item)
        data = var[:]
        if isinstance(data, np.ma.MaskedArray):
            if not allow_masked:
                try:
                    fill_value = var.getncattr("_FillValue")
                except AttributeError:
                    if item == "quality_flags":
                        fill_value = 0
                data = data.filled(fill_value=fill_value)

        return data

    def get_variable_values(self, item, allow_masked=False):
        return self.get_variable_data(
            item, allow_masked=allow_masked
        )  # there is no distinction between data and values in netCDF4

    @classmethod
    def create_dataset_from_parts(
        cls, *, variables, coordinates, attributes, base_ds=None
    ):
        unique_name = f"dataset_{time.time_ns()}_{os.getpid()}.nc"
        ds: netCDF4.Dataset = base_ds or netCDF4.Dataset(
            unique_name, "w", format="NETCDF4", memory=4000 * 4000 * 4
        )
        for name, coordinate in coordinates.items():
            ds.createDimension(name, len(coordinate))
            ds.createVariable(name, coordinate.dtype, (name,))
            ds.variables[name] = coordinate
        for name, variable in variables.items():
            for i, dim in enumerate(variable.dims):
                ds.createDimension(dim, variable.values.shape[i])

            # The "_FillValue" attribute is handled as a special case by netCDF4 and
            # must therefore be set in the constructor of the variable
            fill_value = variable.attrs.pop("_FillValue")
            new_var = ds.createVariable(
                name,
                variable.values.dtype,
                variable.dims,
                fill_value=fill_value,
            )
            for attr_key, attr_val in variable.attrs.items():
                new_var.setncattr(attr_key, attr_val)
            new_var[:] = variable.values

        for attr_key, attr_val in attributes:
            ds[attr_key] = attr_val

        return ds

    def copy_variable(self, variable: str):
        """
        Returns a reference to the original variable.
        """
        return self.get_variable(variable)


class XarrayDatasetAdapter(DatasetAdapter):
    import xarray as xr

    def __init__(self, ds):
        self.ds = ds

    def get_variable(self, item):  # -> da.Array:
        return self.ds[item]

    def get_variable_data(self, item):  # -> da.Array:
        return self.get_variable(item).data

    def get_variable_values(self, item) -> np.ndarray:
        return self.get_variable(item).values

    def create_dataset_from_parts(
        self, *, variables, coordinates, attributes, base_ds=None
    ):
        def create_variable_from_data(
            data: VariableData,
        ) -> self.xr.DataArray:
            return self.xr.DataArray(
                data.values,
                attrs=data.attrs,
                dims=data.dims,
            )

        if base_ds is not None:
            raise ValueError("The 'base_ds' parameter is not supported for xarray.")

        ds = self.xr.Dataset(
            {name: create_variable_from_data(var) for (name, var) in variables.items()},
            coordinates,
            attributes,
        )

        return ds

    def copy_variable(self, variable: str):
        return copy_variable(self.get_variable(variable))


def is_dask_array(arr) -> bool:
    return "dask" in repr(type(arr))


def is_xarray_type(arr) -> bool:
    return "xarray" in repr(type(arr))


def is_netcdf4_type(arr) -> bool:
    return "netCDF4" in repr(type(arr))
