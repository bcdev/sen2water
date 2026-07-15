#  Copyright (c) 2023 EUMETSAT (contract EUM/CO/19/4600002254/JCh)
#
#  Code authored 2023 by Brockmann Consult GmbH under EUMETSAT Contract
#  EUM/CO/19/4600002254/JCh "Sentinel-3 Synergy Cloud Mask Development".

from abc import ABCMeta
from abc import abstractmethod
from typing import Any

# noinspection PyPackageRequirements
import dask.array as da

# noinspection PyPackageRequirements
import numpy as np


class Algorithm(metaclass=ABCMeta):
    """The algorithm interface.

    @author Martin Böttcher, Ralf Quast
    """

    _dtype: np.dtype
    """The result data type. """
    _m: int
    """The number of input array dimensions. """
    _n: int
    """The number of output array dimensions. """

    def __init__(self, dtype: np.dtype, m: int = 2, n: int = 2):
        """Creates a new algorithm instance.
        :param dtype: The result data type.
        :param m: The number of input array dimensions.
        :param n: The number of output array dimensions.
        """
        self._dtype = dtype
        self._m = m
        self._n = n

    @property
    def dtype(self) -> np.dtype:
        """Returns the result data type."""
        return self._dtype

    @property
    @abstractmethod
    def name(self) -> str:
        """Returns the name of the algorithm."""

    @abstractmethod
    def apply_to(self, *inputs: da.Array, **kwargs) -> da.Array:
        """Applies an algorithm to the given inputs and returns the result.
        :param inputs: The input Dask arrays.
        :param kwargs: Any additional keyword arguments of the algorithm.
        :return: A result Dask array.
        """

    @property
    def meta(self) -> np.ndarray:
        """Returns a surrogate instance of the result, which represents
        the result array type and data type."""
        return np.zeros(0, self._dtype)

    @property
    def nan(self):
        """Returns NaN of appropriate result data type."""
        if self.dtype == np.dtype("int8"):
            return np.int8(-127)
        elif self.dtype == np.dtype("int16"):
            return np.int16(-32767)
        elif self.dtype == np.dtype("int32"):
            return np.int32(-2147483647)
        elif self.dtype == np.dtype("int64"):
            return np.int64(-9223372036854775807)
        else:  # unsigned integral and real types are fine with
            return np.full(1, np.nan, dtype=self.dtype).item()


class BlockAlgorithm(Algorithm, metaclass=ABCMeta):
    """Interface class for an algorithm that constitutes a mapping of
    blocks without overlap between blocks.

    @author Martin Böttcher, Ralf Quast
    """

    def __init__(self, dtype: np.dtype, m: int = 2, n: int = 2):
        """Creates a new algorithm instance.
        :param dtype: The result data type.
        :param m: The number of input array dimensions.
        :param n: The number of output array dimensions.
        """
        super().__init__(dtype, m, n)

    @abstractmethod
    def chunks(self, *inputs: da.Array) -> tuple[int, ...] | None:
        """Returns the shape of computed blocks if `compute_blocks` does
        not preserve shape. If the shape is preserved, `None` is returned."""

    @property
    @abstractmethod
    def created_axes(self) -> list[int] | None:
        """Returns the list of dimensions created by `compute_blocks`.
        If no dimensions are created, `None` is returned."""

    @property
    @abstractmethod
    def dropped_axes(self) -> list[int]:
        """Returns the list of dimensions annihilated by `compute_blocks`.
        If no dimensions are annihilated, an empty list is returned."""

    @abstractmethod
    def compute_block(self, *inputs: np.ndarray, **kwargs) -> np.ndarray:
        """Evaluates the algorithm for a single block of data and returns
        the result.
        :param inputs: The input data.
        :param kwargs: Any additional keyword arguments of the computation.
        :return: The result of the computation.
        """

    def compute_block_typed(self, *inputs: np.ndarray, **kwargs) -> np.ndarray:
        """Evaluates the algorithm for a single block of data and converts
        the result into the correct type, if necessary.
        :param inputs: The input data.
        :param kwargs: Any additional keyword arguments of the computation.
        :return: The result of the computation.
        """
        result = self.compute_block(*inputs, **kwargs)
        assert result.ndim == self._n, (
            f"Algorithm {self.name} returned array of invalid dimension "
            f"{result.ndim} != {self._n}"
        )
        return result.astype(self.dtype, copy=False)

    def apply_to(self, *inputs: da.Array, **kwargs) -> da.Array:
        return da.map_blocks(
            self.compute_block_typed,
            *inputs,
            name=self.name,
            dtype=self.dtype,
            chunks=self.chunks(*inputs),
            drop_axis=self.dropped_axes,
            new_axis=self.created_axes,
            meta=self.meta,
            **kwargs,
        )


class OverlapAlgorithm(Algorithm, metaclass=ABCMeta):
    """Interface class for an algorithm that constitutes a mapping of
    blocks with overlap between blocks.

    @author Martin Böttcher, Ralf Quast
    """

    _overlaps: int | tuple[int, ...] | None
    """The number of elements that each block should share with its
    neighbors. If a tuple then this number can be different for each axis,
    otherwise the same number of elements is shared along all axes. """
    _boundary: np.number | str
    """How to handle the boundaries. Values include 'none' and any constant
    numeric value. The default is 'none'. """
    _trim: bool
    """Whether or not to trim `overlaps` elements from each block after
    calling the compute function. """
    _align_arrays: bool
    """Whether or not to align chunks along equally sized dimensions when
    multiple arrays are provided. """
    _allow_rechunk: bool
    """Allows re-chunking, otherwise chunk sizes need to match, and core
    dimensions are to consist only of one chunk. """

    def __init__(
        self,
        dtype: np.dtype,
        m: int = 2,
        n: int = 2,
        overlaps: int | tuple[int, ...] = 1,
        boundary: Any | None = None,
        trim: bool = True,
        align_arrays: bool = True,
        allow_rechunk: bool = False,
    ):
        """Creates a new algorithm instance.
        :param dtype: The result data type.
        :param m: The number of input array dimensions.
        :param n: The number of output array dimensions.
        :param overlaps: The number of elements that each block should
        share with its neighbors. If a tuple then this number can be
        different for each axis, otherwise the same number of elements is
        shared along all axes.
        :param boundary: How to handle the boundaries. Values include
        `None` and any constant numeric value. The default is `None`.
        :param trim: Whether to trim `overlaps` elements from each
        block after calling the compute function. Set this to `False` if your
        compute function already does this for you.
        :param align_arrays: Whether to align chunks along equally
        sized dimensions when multiple arrays are provided. This allows for
        larger chunks in some arrays to be broken into smaller ones that
        match chunk sizes in other arrays such that they are compatible for
        block function mapping. If this is false, then an error will be
        thrown if arrays do not already have the same number of blocks in
        each dimension.
        :param allow_rechunk: Allows re-chunking, otherwise chunk sizes
        need to match, and core dimensions are to consist only of one chunk.
        """
        super().__init__(dtype, m, n)
        self._overlaps = overlaps
        if boundary is None:
            boundary = "none"
        self._boundary = boundary
        self._trim = trim
        self._align_arrays = align_arrays
        self._allow_rechunk = allow_rechunk

    @abstractmethod
    def chunks(self, *inputs: da.Array) -> tuple[int, ...] | None:
        """Returns the shape of computed blocks if `compute_blocks` does
        not preserve shape. If the shape is preserved, `None` is returned."""

    @abstractmethod
    def created_axes(self) -> list[int] | None:
        """Returns the list of dimensions created by `compute_blocks`.
        If no dimensions are created, `None` is returned."""

    @abstractmethod
    def dropped_axes(self) -> list[int]:
        """Returns the list of dimensions annihilated by `compute_blocks`.
        If no dimensions are annihilated, an empty list is returned."""

    @abstractmethod
    def compute_block(self, *inputs: np.ndarray, **kwargs) -> np.ndarray:
        """Evaluates the algorithm for a single block of data and returns
        the result.
        :param inputs: The input data.
        :param kwargs: Any additional keyword arguments of the computation.
        :return: The result of the computation.
        """

    def compute_block_typed(self, *inputs: np.ndarray, **kwargs) -> np.ndarray:
        """Evaluates the algorithm for a single block of data and converts
        the result into the correct type, if necessary.
        :param inputs: The input data.
        :param kwargs: Any additional keyword arguments of the computation.
        :return: The result of the computation.
        """
        result = self.compute_block(*inputs, **kwargs)
        assert result.ndim == self._n, (
            f"Algorithm {self.name} returned array of invalid dimension "
            f"{result.ndim} != {self._n}"
        )
        return result.astype(self.dtype, copy=False)

    def apply_to(self, *inputs: da.Array, **kwargs) -> da.Array:
        return da.map_overlap(
            self.compute_block_typed,
            *inputs,
            depth=self._overlaps,
            boundary=self._boundary,
            trim=self._trim,
            align_arrays=self._align_arrays,
            allow_rechunk=self._allow_rechunk,
            name=self.name,
            dtype=self.dtype,
            chunks=self.chunks(*inputs),
            drop_axis=self.dropped_axes,
            new_axis=self.created_axes,
            meta=self.meta,
            **kwargs,
        )
