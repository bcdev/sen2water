#!/bin/bash
#set -x
set -e

# Linux script to run the Sen2Water processing chain that generates an MSI L2W from an MSI L1C.

s2wdir=$(dirname $(dirname $(realpath $0)))
echo "Sen2Water 0.6.4 ($s2wdir)"

if [ "$(which gpt)" != "$s2wdir/lib/snap/bin/gpt" ]; then
    echo "setting up environment ..."
#    export INSTALL4J_JAVA_HOME=${s2wdir}/lib/jre
#    export INSTALL4J_JAVA_HOME_OVERRIDE=$s2wdir/lib/jre
#    export LD_LIBRARY_PATH="${s2wdir}/lib/snap/snap/modules/lib/amd64:$LD_LIBRARY_PATH"
#    export PYTHONPATH=${s2wdir}/lib/polymer:${s2wdir}/lib/msiresampling:$PYTHONPATH
#    export PATH="${s2wdir}/bin:${s2wdir}/lib/snap/bin:$PATH"
#    export GPT_ADD_CLASSPATH=$s2wdir/lib/snap/snap/modules:$s2wdir/lib/snap/snap/modules/s2tbx-gdal-reader-9.0.2.1cv.jar:$s2wdir/lib/snap/snap/modules/lib-gdal-9.0.6.1cv.jar
#    . ${s2wdir}/lib/conda/bin/activate
#    export ECCODES_DEFINITION_PATH=${s2wdir}/lib/conda/share/eccodes/definitions
    . $s2wdir/bin/mys2w
fi

# parse parameters
input=
c2rccanc=
acoliteanc=
polymeranc=
dem=
withouttgc=false
withdetfoofilter=
chunksize=305
withoutcleanup=
withtoolbox=false
breakpoint="idepix_elevation idepix_slope_aspect_orientation idepix_classification idepix_filter_buffer idepix_mountain_shadow idepix_cloud_statistics idepix_shadow_clustering idepix_cloud_shadow idepix_combine_flags"
outputdir=
while [ "$1" != "" ]; do
    if [ "$1" = "--c2rccanc" ]; then c2rccanc="$2"; shift 2
    elif [ "$1" = "--acoliteanc" ]; then acoliteanc="$2"; shift 2
    elif [ "$1" = "--polymeranc" ]; then polymeranc="$2"; shift 2
    elif [ "$1" = "--dem" ]; then dem="$2"; shift 2
    elif [ "$1" = "--withouttgc" ]; then withouttgc=true; shift 1
    elif [ "$1" = "--withdetfoofilter" ]; then withdetfoofilter="$1"; shift 1
    elif [ "$1" = "--withoutcleanup" ]; then withoutcleanup=true; shift 1
    elif [ "$1" = "--withtoolbox" ]; then withtoolbox=true; shift 1
    elif [ "$1" = "--outputdir" ]; then outputdir="$2"; shift 2
    elif [ "$1" = "--chunksize" ]; then chunksize="$1 $2"; shift 2
    elif [ "$1" = "--breakpoint" ]; then breakpoint="$2"; shift 2
    elif [ "$1" = "--help" ]; then shift 1
    elif [ "${1:0:1}" = "-" ]; then echo unknown parameter $1; exit 1
    elif [ "$input" = "" ]; then input="$1"; shift 1
    else echo "unknown parameter $1"; exit 1
    fi
done

if [ "$input" = "" ]; then
    echo "sen2water.sh <options> <l1cpath>"
    echo "e.g."
    echo "sen2water.sh S2A_MSIL1C_20240104T103431_N0510_R108_T32UME_20240104T123149.SAFE"
    echo "options"
    echo "--c2rccanc embedded | constant"
    echo "--acoliteanc embedded | constant"
    echo "--polymeranc embedded | nasa | constant"
    echo "--dem 'Copernicus 90m Global DEM' | 'Copernicus 30m Global DEM'"
    echo "--withouttgc"
    echo "--withdetfoofilter"
    echo "--withoutcleanup"
    echo "--chunksize 610 | 1830 | 915 | 366 | 305 | 183 | 122 | 61"
    echo "--breakpoint elevation | slope_aspect orientation | idepix_classification | idepix_filter_buffer | idepix_mountain_shadow | idepix_cloud_statistics | idepix_shadow_clustering | idepix_cloud_shadow | idepix_combine_flags"
    exit 1
fi

# S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
if [ ${input:(-1)} = "/" ]; then input=${input:0:(-1)}; fi
if [[ "$input" =~ ".xml" ]]; then input=$(dirname $input); fi
base=$(basename ${input%.SAFE})
granule=${base:39:5}
p=${base:0:3}
y=${base:11:4}
m=${base:15:2}
d=${base:17:2}
H=${base:20:2}
M=${base:22:2}
S=${base:24:2}
resampled=${base}-resampled.nc
destriped=${base}-resampled_TGC.nc
idepix=${base}-idepix.nc
c2rcc=${base}-c2rcc.nc
acolite_output=${p}_MSI_${y}_${m}_${d}_${H}_${M}_${S}_S2R_L2R.nc
acolite=${base}-acolite.nc
cloudmask=${base}-mask.nc
polymer=${base}-polymer.nc
s2w=${base}-s2w.nc

staticmask=${s2wdir}/auxdata/s2w-mask/${granule:0:2}/s2w-mask-${granule}.tif
if [ ! -e $staticmask ]; then
    staticmask=${s2wdir}/auxdata/s2w-global-mask/${granule:0:2}/s2w-globalmask-${granule}.tif
    if [ ! -e $staticmask ]; then
        echo "missing mask ${s2wdir}/auxdata/s2w-global-mask/${granule:0:2}/s2w-globalmask-${granule}.tif"
        echo "assuming ocean"
        staticmask=ocean
    fi
fi

if $withtoolbox; then
    if [ "$outputdir" == "" ]; then
        cd $(dirname $input)
    else
        mkdir -p $outputdir
        cd $outputdir
    fi
else
    if [ "$outputdir" != "" ]; then
        mkdir -p $outputdir
        cd $outputdir
    fi
fi
echo "working directory $(pwd)"

# check whether we modify the software directory
wd=$(realpath .)
if [ "$s2wdir" = "${wd:0:${#s2wdir}}" ]; then
    echo "$(basename $0) must not be started from within software dir $s2wdir."
    echo "cd to some working directory outside the software installation, please."
    exit 1
fi

if ! ln -sf ${s2wdir}/lib/snap/snap/modules/lib/amd64/libenvironment-variables.so . 2> /dev/null; then
    cp ${s2wdir}/lib/snap/snap/modules/lib/amd64/libenvironment-variables.so .
fi

if $withtoolbox; then
    echo 'Progress[%]: 55.0 : Idepix cloud screening ...'
else
    echo "Idepix cloud screening ..."
fi

if [ "$dem" = "" ]; then
    dem="Copernicus 90m Global DEM"
fi
if [ "$chunksize" = "" ]; then
    blocksize=610
else
    blocksize=$(echo $chunksize|awk '{print $2}')
fi
for bp in $breakpoint; do
    idepix=${base}-${bp}.nc 
    # -J-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=8009 
    time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
         -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
         $s2wdir/etc/idepix-graph.xml -Pdem="$dem" -Pbreakpoint="$bp" $destriped -t $idepix -f NetCDF4-BEAM
    echo $idepix
done

if $withtoolbox; then
    rm -f /tmp/s2woutput/*
    mkdir -p /tmp/s2woutput
    ln -s $(pwd)/$newname /tmp/s2woutput/$newname
fi

if [ "$withoutcleanup" != "true" ]; then
    #rm $resampled $idepix $c2rcc $acolite ${acolite/L2R/L1R} ${polymer/polymer/mask} $polymer
    rm -f ${base}-TGC-parameters.json $destriped 
    rm -f acolite_run*txt ${base}-acolite.log acolite.parameters polymer.parameters
    rm -f libenvironment-variables.so
fi

if $withtoolbox; then
    echo Progress[%%]: 100.0 : done
    echo $(pwd)/$newname
else
    echo $newname
fi
echo "done"
