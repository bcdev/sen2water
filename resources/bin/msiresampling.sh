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
if [ "$(which gpt)" != "$s2wdir/lib/snap/bin/gpt" ]; then
    echo "Environment not set up. Please source mys2w once with"
    echo ". $s2wdir/bin/mys2w"
    exit 1
fi
# check whether we have an input parameter
if [ "$1" = "" ]; then
    echo "msiresampling.sh <l1cpath>"
    echo "e.g."
    echo "msiresampling.sh S2A_MSIL1C_20240104T103431_N0510_R108_T32UME_20240104T123149.SAFE"
    exit 1
fi

# S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
input=$1
base=$(basename ${input%.SAFE})
resampled=${base}-resampled.nc

echo "resampling to 60m ..."

time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py --progress $input $resampled

echo $resampled
echo "done"
