import xarray as xr
import numpy as np
import dask.array as da
import dask
from abc import ABC, abstractmethod

######################################## data ########################################

def source_dt() -> xr.DataTree:
    dims = "y", "x"
    data_array = da.from_array(np.arange(16).reshape(4, 4), chunks=(2,2))
    data = xr.DataArray(data_array, dims=dims)
    return xr.DataTree.from_dict(
        { 
            "/data": data,
        }
    )

######################################## algorithms ########################################

@dask.delayed
def compute_mean(arr):
    return arr.mean()


######################################## processing units ########################################

class ProcessingUnit(ABC):
    @abstractmethod
    # *dt
    def apply(self, *dt: xr.DataTree, **kwargs) -> xr.DataTree: ...

class BlockProcessingUnit(ABC):
    @abstractmethod
    def func(self, *arr: np.ndarray, **kwargs) -> np.ndarray: ...

    def _map_dask_blocks(self, *inputs, **kwargs) -> da.Array:
        return da.map_blocks(self.func, *inputs, **kwargs)

    def apply_internal(self, *inputs: xr.Dataset, **kwargs) -> xr.Dataset:
        return xr.apply_ufunc(self._map_dask_blocks, *inputs, dask="allowed", **kwargs)

    #def apply(dt: xr.DataTree, **kwargs) -> xr.DataTree:
    #    # Alternative: extract and reconstruct data tree
    #    xr.apply_ufunc(self._map_dask_blocks, dt, kwargs=kwargs) # TODO likely must be a data set/array


class ReplaceWithMaxPU(BlockProcessingUnit):
    def replace_with_max(self, arr: np.ndarray, **kwargs) -> np.ndarray:
        new_arr = np.empty(arr.shape, dtype=arr.dtype)
        new_arr.fill(arr.max())
        return new_arr

    func = replace_with_max

    def apply(self, dt: xr.DataTree, **kwargs) -> xr.DataTree:
        input = dt["/data"]
        output = self.apply_internal(input, kwargs={**kwargs, "dtype": float})
        tree = xr.DataTree.from_dict({"/max": output})
        return tree


class GlobalMeanPU(BlockProcessingUnit):
    def block_mean(self, arr: np.ndarray, **kwargs) -> np.ndarray:
        single_element_array = np.empty(np.ones(arr.shape), dtype=arr.dtype)
        single_element_array.fill(arr.mean())
        return single_element_array

    func = block_mean

    def global_mean(self, arr: da.Array, **kwargs) -> float:
        return arr.values.mean()

    def apply(self, dt: xr.DataTree, **kwargs):
        means = self.apply_internal(self.block_mean, dt["/max"].data, kwargs={**kwargs, "dtype": float})
        global_mean_data = dask.delayed(self.global_mean)(means)
        mean = da.from_delayed(global_mean_data, shape=(1, ), dtype=float)
        # TODO check this is a proper data tree that can be stored
        tree = xr.DataTree.from_dict({"/mean": xr.DataArray(mean)})
        return tree

class SubtractScalarPU(BlockProcessingUnit):
    def subtract(self, arr, scalar=None, **kwargs) -> np.ndarray:
        return arr - scalar

    func = subtract

    def apply(self, dt: xr.DataTree, dt_mean: xr.DataTree, **kwargs) -> xr.DataTree:
        input = dt["/max"]
        global_mean = dt_mean["/mean"].data  # TODO we have to unwrap keyword args ourselves here
        output = self.apply_internal(input, kwargs={**kwargs, "scalar":global_mean, "dtype": float})
        tree = xr.DataTree.from_dict({ "/sub": xr.DataArray(output) })
        return tree

######################################## structure ########################################


def main():
    dt = source_dt()
    dt_replaced = ReplaceWithMaxPU().apply(dt)
    print(f"{dt_replaced=}")
    dt_mean = GlobalMeanPU().apply(dt_replaced)
    print(f"{dt_mean=}")
    dt_subtracted = SubtractScalarPU().apply(dt_replaced, dt_mean)
    print(f"{dt_subtracted=}")
    dt_subtracted["/sub"].data.visualize("unoptimized.svg")
    #dt_subtracted["/sub"].data.visualize("unoptimized-cyto.svg", engine="cytoscape")
    dt_subtracted["/sub"].data.visualize("optimized.svg", optimize_graph=True)
    #dt_subtracted["/sub"].data.visualize("optimized-cyto.svg", optimize_graph=True, engine="cytoscape")
    #print(f"{dt_subtracted['/sub'].compute()=}")
    print("Done")

if __name__ == "__main__":
    main()
