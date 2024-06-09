#!/bin/bash
#set -x
set -e

s2wdir=$(dirname $(dirname $(realpath $0)))

# check whether we modify the software directory
wd=$(realpath .)
if [ "$s2wdir" = "${wd:0:${#s2wdir}}" ]; then
    echo "$(basename $0) must not be started from within software dir $s2wdir."
    echo "cd to some working directory outside the software installation, please."
    exit 1
fi
if [ "$(which python)" != "$s2wdir/lib/conda/bin/python" ]; then
    echo "Environment not set up. Please source mys2w once with"
    echo ". $s2wdir/bin/mys2w"
    exit 1
fi
# check whether we have an input parameter
if [ "$1" = "" ]; then
    echo "hrocresampling.sh <l1cpath> [<hrocmask>]"
    echo "e.g."
    echo "hrocresampling.sh S2A_MSIL1C_20240104T103431_N0510_R108_T32UME_20240104T123149.SAFE raster_mask_tile_32UME_4nws_32632.tif"
    exit 1
fi

# S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
input=$1
base=$(basename ${input%.SAFE})
granule=${base:39:5}
if [ "$2" != "" ]; then
    # parameter
    hrocmask=$2
elif [ -e ${s2wdir}/auxdata/s2w-mask/${granule:0:2}/s2w-mask-${granule}.tif ]; then
    # European coastal ocean/inlandwater/land mask
    hrocmask=${s2wdir}/auxdata/s2w-mask/${granule:0:2}/s2w-mask-${granule}.tif
else
    # Global ocean/inlandwater/land mask
    hrocmask=${s2wdir}/auxdata/s2w-global-mask/${granule:0:2}/s2w-globalmask-${granule}.tif
fi
resampled=${base}-hrocresampled.nc

echo "resampling to 60m ..."

time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py \
     --ancillary msl --ancillary tco3 --ancillary tcwv --ancillary u10 --ancillary v10 \
     --hrocmask $hrocmask --progress $input $resampled

echo $resampled
echo "done"
