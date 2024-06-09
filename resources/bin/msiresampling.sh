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
# check whether we have a proper env
if [ "$(which python)" != "$s2wdir/lib/conda/bin/python" ]; then
    echo "Environment not set up. Please source mys2w once with"
    echo ". $s2wdir/bin/mys2w"
    exit 1
fi

time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py $@

echo "done"
