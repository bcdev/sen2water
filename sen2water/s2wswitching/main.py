# -*- coding: utf-8 -*-

"""Command line client of Sen2Water's switching processor that combines AC results into L2W"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "TBD"
__version__ = "0.51"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

# changes in 0.51:
# separate to_netcdf from computation of statistics, attempt to avoid infinite run of compute

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
@click.option("--chunksize",
              type=click.Choice(['1830', '915', '610', '366', '305', '183', '122', '61']),
              default='610')
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
    chunksize: str,
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
                chunksize,
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
    chunksize: str,
    with_copyinputs: bool,
) -> int:
    if chunksize:
        blocksize = int(chunksize)
    else:
        blocksize = 610
    """Converts paths to xarray Datasets, writes output Dataset to file"""
    try:
        logger.info("opening inputs")
        input_id = os.path.basename(resampled).replace(".zip", "").replace(".SAFE", "")
        resampled_ds = xr.open_dataset(
            resampled, chunks={"y": blocksize, "x": blocksize}, mask_and_scale=False
        )
        idepix_ds = xr.open_dataset(
            idepix, chunks={"y": blocksize, "x": blocksize}, mask_and_scale=False
        )
        c2rcc_ds = xr.open_dataset(
            c2rcc, chunks={"y": blocksize, "x": blocksize}, mask_and_scale=False
        )
        acolite_ds = xr.open_dataset(
            acolite, chunks={"y": blocksize, "x": blocksize}, mask_and_scale=False
        )
        polymer_ds = xr.open_dataset(
            polymer, chunks={"height": blocksize, "width": blocksize}, mask_and_scale=True  # Polymer does not use NaN as fill value
        )
        s2wmask_ds = rio.open_rasterio(
            s2wmask, chunks={"y": blocksize, "x": blocksize}, mask_and_scale=False
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
#        write_to_netcdf = output_ds.to_netcdf(
        output_ds.to_netcdf(
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
#            compute=False,
        )
        logger.info("data written")
        logger.info("computing statistics")
        count_pixels = S2wStatistics.count_pixels(
            output_ds["pixel_class"].data, s2wmask_ds["s2wmask"].data
        )
#        logger.info("starting computation")
#        statistics = dask.compute([*count_pixels, write_to_netcdf])
        statistics = dask.compute(*count_pixels)
        logger.info("adding statistics")
#        statistics_attrs = dict(zip(S2wStatistics.attributes(), statistics[0][:-1]))
        statistics_attrs = dict(zip(S2wStatistics.attributes(), statistics))
        statistics_ds = xr.Dataset(attrs=statistics_attrs)
        statistics_ds.to_netcdf(output, mode="a")
        logger.info(f"output {output} written")
        return 0
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return 1


if __name__ == "__main__":

    import code, traceback, signal
    def debug(sig, frame):
        """Interrupt running process, and provide a python prompt for
        interactive debugging."""
        d={'_frame':frame}         # Allow access to frame object.
        d.update(frame.f_globals)  # Unless shadowed by global
        d.update(frame.f_locals)
        i = code.InteractiveConsole(d)
        message  = "Signal received : entering python shell.\nTraceback:\n"
        message += ''.join(traceback.format_stack(frame))
        i.interact(message)
    def listen():
        signal.signal(signal.SIGQUIT, debug)  # Register handler
    listen()

    sys.exit(run())
