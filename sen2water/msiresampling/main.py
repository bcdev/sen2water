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

"""
Open points
- find out how to split a stacked map_blocks output to avoid repeated detector majority computation
- check chunking parameter in some calls
"""

import sys
import click
import traceback
import os.path
from typing import List
import dask.array
import xarray as xr
from sen2water.eoutils.eologging import logger
from sen2water.eoutils.eoscheduler import Scheduler
from sen2water.eoutils.eoprofiling import Profiling
from sen2water.eoutils.eoprogress import Progress
from sen2water.msiresampling.resamplingoperator import ResamplingOperator


@click.command()
@click.argument("l1c")
@click.argument("output")
@click.option("--resolution",
              type=click.Choice(['10', '20', '60']),
              default='60')
@click.option("--downsampling",
              type=click.Choice(['detectormean', 'first', 'min', 'max', 'mean', 'median']),
              default='detectormean')
@click.option("--flagdownsampling",
              type=click.Choice(['first', 'flagand', 'flagor', 'flagmedianand', 'flagmedianor']),
              default='first')
@click.option("--upsampling",
              type=click.Choice(['nearest', 'bilinear', 'bicubicnotyetsupported']),
              default='bilinear')
@click.option("--ancillary",
              type=click.Choice(['msl', 'tco3', 'tcwv', 'u10', 'v10', 'aod550']),
              multiple=True)
@click.option("--withmasterdetfoo", is_flag=True)
@click.option("--scheduler",
              type=click.Choice(["synchronous", "threads", "processes"]),
              default="threads")
@click.option("--profiling")
@click.option("--progress", is_flag=True)
def run(
    l1c: str,
    output: str,
    resolution: str,
    downsampling: str,
    flagdownsampling: str,
    upsampling: str,
    ancillary: List[str],
    withmasterdetfoo: bool,
    scheduler: str,
    profiling: str,
    progress: bool,
) -> int:
    """Selects scheduler and optionally uses profiling"""
    if profiling and scheduler != "synchronous":
        print("profiling uses synchronous scheduler")
        scheduler = "synchronous"
    with Scheduler(scheduler):
        with Profiling(profiling):
            code = _run(
                l1c,
                output,
                resolution,
                downsampling,
                flagdownsampling,
                upsampling,
                ancillary,
                withmasterdetfoo,
                progress,
            )
    return code


def _run(
    l1c: str,
    output: str,
    resolution: str,
    downsampling: str,
    flagdownsampling: str,
    upsampling: str,
    ancillary: List[str],
    withmasterdetfoo: bool,
    progress: bool,

) -> int:
    """Converts paths to xarray Datasets, writes output Dataset to file"""
    try:
        resampler = ResamplingOperator()
        logger.info("opening inputs")
        input_id = os.path.basename(l1c).replace(".zip", "").replace(".SAFE", "")
        l1c_ds = xr.open_dataset(l1c, chunks=resampler.preferred_chunks(), engine="safe_msi_l1c", merge_flags=True)
        logger.info("starting computation")
        output_ds = resampler.run(
            l1c_ds, int(resolution), downsampling, flagdownsampling, upsampling, ancillary, withmasterdetfoo,
            merge_flags=True
        )
        with Progress(progress):
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
    sys.exit(run())
