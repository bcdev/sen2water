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

def replace_with_max(arr: np.ndarray) -> np.ndarray: # unchanged shape
    new_arr = arr.copy()
    new_arr[:] = arr.max()
    return new_arr

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
    
    #def apply(dt: xr.DataTree, **kwargs) -> xr.DataTree:
    #    # Alternative: extract and reconstruct data tree
    #    xr.apply_ufunc(self._map_dask_blocks, dt, kwargs=kwargs) # TODO likely must be a data set/array


class ReplaceWithMaxPU(BlockProcessingUnit):
    def replace_with_max(self, arr, **kwargs) -> np.ndarray:
        return replace_with_max(arr)

    func = replace_with_max

    def apply(self, dt: xr.DataTree, **kwargs) -> xr.DataTree:
        transformed_data = xr.apply_ufunc(self._map_dask_blocks, dt["/data"], dask="allowed", kwargs=kwargs)
        new_dt = xr.DataTree.from_dict(
            {
                #"/data": xr.Dataset({"data": transformed_data})
                "/max": transformed_data,
            }
        )
        return new_dt

class GlobalMeanPU(ProcessingUnit):
    def apply(self, dt: xr.DataTree, **kwargs):
        new_dt = dt.copy()
        computed_mean = compute_mean(dt["/max"].data)
        new_dt["/mean"] = da.from_delayed(computed_mean, shape=(1, ), dtype=float)
        return new_dt

class SubtractScalarPU(BlockProcessingUnit):
    def subtract(self, arr, **kwargs) -> np.ndarray:
        return arr - kwargs["scalar"]

    func = subtract

    def apply(self, dt: xr.DataTree, dt_mean: xr.DataTree, **kwargs) -> xr.DataTree:
        #transformed_data = xr.apply_ufunc(self._map_dask_blocks, dt["/max"], dask="allowed", kwargs=kwargs | {"scalar": dt_mean["/mean"].data, "dtype": float})
        transformed_data = self._map_dask_blocks(dt["/max"].data, scalar=dt_mean["/mean"].data, dtype=float, **kwargs)
        new_dt = xr.DataTree.from_dict(
            {
                "/sub": xr.DataArray(transformed_data)
            }
        )
        return new_dt

######################################## structure ########################################


def main():
    dt = source_dt()
    dt_replaced = ReplaceWithMaxPU().apply(dt)
    dt_mean = GlobalMeanPU().apply(dt_replaced)
    print(f"{dt_mean=}")
    dt_subtracted = SubtractScalarPU().apply(dt_replaced, dt_mean)
    print(f"{dt_subtracted=}")
    dt_subtracted["/sub"].data.visualize("unoptimized.svg")
    dt_subtracted["/sub"].data.visualize("unoptimized-cyto.svg", engine="cytoscape")
    dt_subtracted["/sub"].data.visualize("optimized.svg", optimize_graph=True)
    dt_subtracted["/sub"].data.visualize("optimized-cyto.svg", optimize_graph=True, engine="cytoscape")
    print(f"{dt_subtracted['/sub'].compute()=}")
    print("Done")

if __name__ == "__main__":
    main()
