# -*- coding: utf-8 -*-

"""Command line client of Sen2Water's switching processor that combines AC results into L2W"""

__author__ = "Martin Böttcher, Brockmann Consult GmbH"
__copyright__ = "Copyright 2023, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Production"

# changes in 0.51:
# separate to_netcdf from computation of statistics, attempt to avoid infinite run of compute
# changes in 0.52
# use h5netcdf engine for writing to avoid deadlock with reading in netcdf4
# changes in 0.53
# add dummy band to acolite in case a band is missing
# changes in 0.6
# license header

import sys
import traceback
import click
import warnings
import os.path
import dask
import xarray as xr
import rioxarray as rio
from sen2water.eoutils.eologging import logger
from sen2water.eoutils.eoscheduler import Scheduler
from sen2water.eoutils.eoprofiling import Profiling
from sen2water.s2wswitching.switchingoperator import SwitchingProcessor
from sen2water.s2wswitching.statistics import S2wStatistics
from sen2water.s2wswitching.acolitebands import AcoliteDataset


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
    sys.exit(code)


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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
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
            acolite_ds = AcoliteDataset(acolite_ds)
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
        try:
            os.unlink(output)  # avoids deadlock in netcdf4 lib
        except:
            pass
        logger.info("starting computation")
        output_ds.to_netcdf(
            output,
            engine="h5netcdf",
            encoding={
                **{
                    var: {
                        "zlib": True,
                        "complevel": 5,
                        "chunksizes": output_ds[var].data.chunksize,
                    }
                    for var in output_ds.data_vars
                }
            }
        )
        logger.info("computing statistics")
        count_pixels = S2wStatistics.count_pixels(
            output_ds["pixel_class"].data, s2wmask_ds["s2wmask"].data
        )
        statistics = dask.compute(*count_pixels)
        logger.info("adding statistics")
        statistics_attrs = dict(zip(S2wStatistics.attributes(), statistics))
        statistics_ds = xr.Dataset(attrs=statistics_attrs)
        statistics_ds.to_netcdf(output, mode="a")
        logger.info(f"output {output} written")
        return 0
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        return 1


if __name__ == "__main__":

    try:  # only valid in Linux
        import threading, signal
        def thread_dump(_sig, _frame):
            for th in threading.enumerate():
                print(th)
                traceback.print_stack(sys._current_frames()[th.ident])
                print()
        def listen():
            signal.signal(signal.SIGQUIT, thread_dump)
        listen()
    except:
        pass
    
    run()
