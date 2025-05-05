# -*- coding: utf-8 -*-

"""
Command line client of the msiresampling processor.
msiresampling obeys detectors when interpolating viewing angles
and interpolates atmospheric auxiliary data using the correct coarse grid.
It uses dask threads for parallelisation and memory management.
"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 1.1:
# ...

import sys
import warnings
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
from sen2water.msiresampling.hrocmask import HrocMask
from sen2water.msiresampling.resamplingoperator import ResamplingOperator

@click.command()
@click.argument("l1c")
@click.argument("output", required=False)
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
@click.option("--withdetfoofilter", is_flag=True)
@click.option("--chunksize",
              type=click.Choice(['1830', '915', '610', '366', '305', '183', '122', '61']),
              default='610')
@click.option("--hrocmask")
@click.option("--scheduler",
              type=click.Choice(["synchronous", "threads", "processes"]),
              default="threads")
@click.option("--profiling")
@click.option("--progress/--noprogress", is_flag=True, default=True)
def run(
    l1c: str,
    output: str,
    resolution: str,
    downsampling: str,
    flagdownsampling: str,
    upsampling: str,
    ancillary: List[str],
    withdetfoofilter: bool,
    chunksize: str,
    hrocmask: str,
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
                withdetfoofilter,
                chunksize,
                hrocmask,
                progress,
            )
    sys.exit(code)


def _run(
    l1c: str,
    output: str,
    resolution: str,
    downsampling: str,
    flagdownsampling: str,
    upsampling: str,
    ancillary: List[str],
    withdetfoofilter: bool,
    chunksize: str,
    hrocmask: str,
    progress: bool,

) -> int:
    """Converts paths to xarray Datasets, writes output Dataset to file"""
    try:
        logger.info("starting resampling")
        resampler = ResamplingOperator(int(chunksize) * 60)
        input_id = os.path.basename(l1c).replace(".zip", "").replace(".SAFE", "")
        if not output:
            output = f"{input_id}-resampled.nc"
        logger.info("opening inputs")
        l1c_ds = xr.open_dataset(
            l1c, chunks=resampler.preferred_chunks(), engine="safe_msi_l1c", merge_flags=True
        )
        # open static mask before resampling to fail early in case it is missing
        if hrocmask and hrocmask != "ocean":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                hroc_mask_ds = xr.open_dataset(
                    hrocmask, chunks={ "y": 610, "x": 610 }, engine="rasterio"
                )
        else:
            hroc_mask_ds = None
            HROC_MASK_OCEAN = 2
        # change default of ancillary interpolation in case of hroc, using static mask as marker
        if hrocmask and not ancillary:
            ancillary = ['msl', 'tco3', 'tcwv', 'u10', 'v10']
        logger.info("preparing computation")
        output_ds = resampler.run(
            l1c_ds,
            int(resolution),
            downsampling,
            flagdownsampling,
            upsampling,
            ancillary,
            withdetfoofilter,
            merge_flags=True
        )
        if hrocmask:
            if hroc_mask_ds:
                output_ds = HrocMask().run(output_ds, hroc_mask_ds)
            else:
                output_ds = HrocMask().run_with_constant(output_ds, HROC_MASK_OCEAN)
        logger.info("starting computation")
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
    run()
