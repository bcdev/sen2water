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
import os.path
import dask
import xarray as xr
import rioxarray as rio
from sen2water.eoutils.eologging import logger
from sen2water.eoutils.eoscheduler import Scheduler
from sen2water.eoutils.eoprofiling import Profiling
from sen2water.s2wswitching.switchingoperator import SwitchingProcessor
from sen2water.s2wswitching.statistics import S2wStatistics


@click.command()
@click.argument("resampled")
@click.argument("idepix")
@click.argument("c2rcc")
@click.argument("acolite")
@click.argument("polymer")
@click.argument("s2wmask")
@click.argument("output")
@click.option("--copyinputs", "with_copyinputs", is_flag=True, default=False)
@click.option(
    "--scheduler",
    type=click.Choice(["synchronous", "threads", "processes"]),
    default="threads",
)
@click.option("--profiling")
def run(
    resampled: str,
    idepix: str,
    c2rcc: str,
    acolite: str,
    polymer: str,
    s2wmask: str,
    output: str,
    with_copyinputs: bool,
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
                resampled,
                idepix,
                c2rcc,
                acolite,
                polymer,
                s2wmask,
                output,
                with_copyinputs,
            )
    return code


def _run(
    resampled: str,
    idepix: str,
    c2rcc: str,
    acolite: str,
    polymer: str,
    s2wmask: str,
    output: str,
    with_copyinputs: bool,
) -> int:
    """Converts paths to xarray Datasets, writes output Dataset to file"""
    try:
        logger.info("opening inputs")
        input_id = os.path.basename(resampled).replace(".zip", "").replace(".SAFE", "")
        resampled_ds = xr.open_dataset(
            resampled, chunks={"y": 610, "x": 610}, mask_and_scale=False
        )
        idepix_ds = xr.open_dataset(
            idepix, chunks={"y": 610, "x": 610}, mask_and_scale=False
        )
        c2rcc_ds = xr.open_dataset(
            c2rcc, chunks={"y": 610, "x": 610}, mask_and_scale=False
        )
        acolite_ds = xr.open_dataset(
            acolite, chunks={"y": 610, "x": 610}, mask_and_scale=False
        )
        polymer_ds = xr.open_dataset(
            polymer, chunks={"height": 610, "width": 610}, mask_and_scale=True  # Polymer does not use NaN as fill value
        )
        s2wmask_ds = rio.open_rasterio(
            s2wmask, chunks={"y": 610, "x": 610}, mask_and_scale=False
        ).to_dataset(name="s2wmask")
        logger.info("inputs opened")
        output_ds = SwitchingProcessor().run(
            resampled=resampled_ds,
            idepix=idepix_ds,
            c2rcc=c2rcc_ds,
            acolite=acolite_ds,
            polymer=polymer_ds,
            s2wmask=s2wmask_ds,
            input_id=input_id,
            with_copyinputs=with_copyinputs,
        )
        write_to_netcdf = output_ds.to_netcdf(
            output,
            encoding={
                **{
                    var: {
                        "zlib": True,
                        "complevel": 5,
                        "chunksizes": output_ds[var].data.chunksize,
                    }
                    for var in output_ds.data_vars
                }
            },
            compute=False,
        )
        logger.info("writing prepared")
        count_pixels = S2wStatistics.count_pixels(
            output_ds["pixel_class"].data, s2wmask_ds["s2wmask"].data
        )
        logger.info("starting computation")
        statistics = dask.compute([write_to_netcdf, *count_pixels])
        logger.info("adding statistics")
        statistics_attrs = dict(zip(S2wStatistics.attributes(), statistics[0][1:]))
        statistics_ds = xr.Dataset(attrs=statistics_attrs)
        statistics_ds.to_netcdf(output, mode="a")
        logger.info(f"output {output} written")
        return 0
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return 1


if __name__ == "__main__":
    sys.exit(run())
