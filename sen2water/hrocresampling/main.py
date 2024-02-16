# -*- coding: utf-8 -*-

"""Command line client of Sen2Water's switching processor that combines AC results into L2W"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.5"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 1.1:
# ...

import sys
import click
import traceback

import dask.array
import xarray as xr
from sen2water.eoutils.eologging import logger
from sen2water.eoutils.eoscheduler import Scheduler
from sen2water.eoutils.eoprofiling import Profiling
from sen2water.hrocresampling.hrocmask import HrocMask
from sen2water.msiresampling.resamplingoperator import ResamplingOperator


@click.command()
@click.argument("l1c")
@click.argument("hrocmask")
@click.argument("output")
@click.option(
    "--scheduler",
    type=click.Choice(["synchronous", "threads", "processes"]),
    default="threads",
)
@click.option("--profiling")
def run(
    l1c: str,
    hrocmask: str,
    output: str,
    scheduler: str,
    profiling: str,
) -> int:
    """Selects scheduler and optionally uses profiling"""
    if profiling and scheduler != "synchronous":
        print("profiling uses synchronous scheduler")
        scheduler = "synchronous"
    with Scheduler(scheduler):
        with Profiling(profiling):
            code = _run(
                l1c,
                hrocmask,
                output,
            )
    return code


def _run(
    l1c: str,
    hrocmask: str,
    output: str,
) -> int:
    """Converts paths to xarray Datasets, writes output Dataset to file"""
    try:
        logger.info("opening inputs")
        resampler = ResamplingOperator()
        l1c_ds = xr.open_dataset(l1c, chunks=resampler.preferred_chunks(), engine="safe_msi_l1c", merge_flags=True)
        hroc_mask_ds = xr.open_dataset(hrocmask, chunks={ "y": 610, "x": 610 }, engine="rasterio")

        logger.info("starting computation")
        intermediate_ds = resampler.run(l1c_ds, 60, "mean", "first", "bilinear", merge_flags=True)
        output_ds = HrocMask().run(intermediate_ds, hroc_mask_ds)

        output_ds.to_netcdf(
            output,
            encoding={
                **{
                    var: {
                        "zlib": True,
                        "complevel": 5,
                        "chunksizes": output_ds[var].data.chunksize,
                    }
                    for var in output_ds.data_vars if isinstance(output_ds[var].data, dask.array.Array)
                }
            }
        )
        logger.info(f"output {output} written")
        return 0
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return 1


if __name__ == "__main__":
    run()
