# -*- coding: utf-8 -*-

"""Generic operator and algorithm interfaces and abstract classes"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

from typing import Any
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as nt
import dask.array as da
import xarray as xr
from sen2water.eoutils.eologging import logger


class Operator(ABC):
    log = logger

    @abstractmethod
    def run(self, *inputs: xr.Dataset, **kwargs: Any) -> xr.Dataset: ...


class Algorithm(ABC):
    @abstractmethod
    # def process_image(self, *inputs: da.Array, **kwargs: Any) -> da.Array:
    def apply(self, *inputs: da.Array, **kwargs: Any) -> da.Array: ...


class BlockAlgorithm(Algorithm):
    # def process_image(self, *inputs: da.Array, dtype: nt.DTypeLike = float, **kwargs: Any) -> da.Array:
    def apply(
        self, *inputs: da.Array, dtype: nt.DTypeLike = float, **kwargs: Any
    ) -> da.Array:
        return da.map_blocks(
            self.func, *inputs, dtype=dtype, meta=np.array((), dtype=dtype), **kwargs
        )

    @abstractmethod
    # def process_block(self, *inputs: np.ndarray, **kwargs: Any) -> np.ndarray:
    def func(self, *inputs: np.ndarray, **kwargs: Any) -> np.ndarray: ...


class OverlapAlgorithm(Algorithm):
    # def process_image(self, *inputs: da.Array, dtype: nt.DTypeLike = float, depth: int = 1, **kwargs: Any) -> da.Array:
    def apply(
        self,
        *inputs: da.Array,
        dtype: nt.DTypeLike = float,
        depth: int = 1,
        **kwargs: Any
    ) -> da.Array:
        return da.map_overlap(
            self.func,
            *inputs,
            depth=depth,
            dtype=dtype,
            meta=np.array((), dtype=dtype),
            **kwargs
        )

    @abstractmethod
    # def process_block(self, *inputs: np.ndarray, **kwargs: Any) -> np.ndarray:
    def func(self, *inputs: np.ndarray, **kwargs: Any) -> np.ndarray: ...
