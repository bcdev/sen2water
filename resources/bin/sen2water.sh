#!/bin/bash
#set -x
set -e

s2wdir=$(dirname $(dirname $(realpath $0)))
echo "Sen2Water 0.6.1 ($s2wdir)"

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
chunksize=
withoutcleanup=
withtoolbox=false
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
acolite=${p}_MSI_${y}_${m}_${d}_${H}_${M}_${S}_S2R_L2R.nc
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
    echo 'Progress[%]: 5.0 : resampling to 60m ...'
else
    echo "resampling to 60m ..."
fi

time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py \
     --ancillary msl --ancillary tco3 --ancillary tcwv --ancillary u10 --ancillary v10 \
     $withdetfoofilter $chunksize $input $resampled

echo $resampled
echo

if $withouttgc; then
    destriped=${resampled}
else
    if $withtoolbox; then
        echo 'Progress[%]: 25.0 : TOA glint correction ...'
    else
        echo "TOA glint correction"
    fi

    time python $s2wdir/lib/acolite/tgcparameters.py ${input}
    time python -u $s2wdir/lib/acolite/tgc2.py --acolite_path $s2wdir/lib/acolite --input ${resampled} --output . \
           --glint_threshold 0.0 --scaling True --aot_min 0.1 --process_high_sza False --estimate False \
           --grid_write True --grid_files ${base}-TGC-parameters.json --verbosity 2
    if [ -e ${destriped} ]; then
        echo $destriped
    else
        echo "*** no ${destriped} found, using resampled instead"
        destriped=${resampled}
    fi
    echo
fi

if $withtoolbox; then
    echo 'Progress[%]: 45.0: ACOLITE atmospheric correction ...'
else
    echo "ACOLITE atmospheric correction ..."
fi

if [ "$acoliteanc" == "constant" ]; then
    s2auxiliarydefault=False
else
    s2auxiliarydefault=True
fi

cat $s2wdir/etc/acolite.parameters | sed \
    -e s,S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc,$destriped, \
    -e s,s2_auxiliary_default=True,s2_auxiliary_default=${s2auxiliarydefault}, \
    > acolite.parameters
time python -u $s2wdir/lib/acolite/launch_acolite.py --nogfx --cli --settings=acolite.parameters > ${base}-acolite.log
if [ ! -e $acolite ]; then
    echo "*** ACOLITE failed, retrying without TGC ..."
    cat $s2wdir/etc/acolite.parameters | sed \
        -e s,S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc,$resampled, \
        -e s,s2_auxiliary_default=True,s2_auxiliary_default=${s2auxiliarydefault}, \
        > acolite.parameters
    time python -u $s2wdir/lib/acolite/launch_acolite.py --nogfx --cli --settings=acolite.parameters > ${base}-acolite.log
    if [ ! -e $acolite ]; then
        echo "*** ACOLITE failed"
        exit 1
    fi
    destriped=$resampled
fi

echo $acolite
echo

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
time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
    -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
    $s2wdir/etc/idepix-graph.xml -Pdem="$dem" $destriped -t $idepix -f NetCDF4-BEAM

echo $idepix
echo

if $withtoolbox; then
    echo 'Progress[%]: 65.0 : C2RCC atmospheric correction ...'
else
    echo "C2RCC atmospheric correction ..."
fi

if [ "$c2rccanc" == "embedded" ]; then
    useEcmwfAuxData=true
else
    useEcmwfAuxData=false
fi
time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
    -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
    $s2wdir/etc/c2rcc-graph.xml -Pdem="$dem" -PuseEcmwfAuxData="$useEcmwfAuxData" $destriped -t $c2rcc -f NetCDF4-BEAM

echo $c2rcc
echo

if $withtoolbox; then
    echo 'Progress[%]: 75.0 : POLYMER atmospheric correction ...'
else
    echo "POLYMER atmospheric correction ..."
fi

time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
    -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
    $s2wdir/etc/polymer-mask-graph.xml $idepix -t $cloudmask -f NetCDF4-BEAM

if [ "$polymeranc" == "" ]; then
    polymeranc=embedded
fi

cat $s2wdir/etc/polymer.parameters | sed -e s/S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_mask.nc/$cloudmask/ > polymer.parameters
rm -f $polymer
time python -u $s2wdir/lib/polymer/run-polymer.py polymer.parameters "$polymeranc" "$s2wdir/auxdata/dem/$dem" $resampled $polymer

echo $polymer
echo

if $withtoolbox; then
    echo 'Progress[%]: 85.0 : Sen2Water switching and output formatting ...'
else
    echo "Sen2Water switching and output formatting ..."
fi

time python -u $s2wdir/lib/msiresampling/sen2water/s2wswitching/main.py \
     $chunksize \
     $destriped $idepix $c2rcc $acolite $polymer $staticmask $s2w

if $withtoolbox; then
    echo 'Progress[%]: 95.0 : renaming output ...'
else
    echo "renaming output ..."
fi

newname=$(ncdump -h $s2w | grep ':id' | cut -d '"' -f 2)
mv $s2w $newname

if $withtoolbox; then
    rm -f /tmp/s2woutput/*
    mkdir -p /tmp/s2woutput
    ln -s $(pwd)/$newname /tmp/s2woutput/$newname
fi

if [ "$withoutcleanup" != "true" ]; then
    rm $resampled $idepix $c2rcc $acolite ${acolite/L2R/L1R} ${polymer/polymer/mask} $polymer
    rm -f ${base}-TGC-parameters.json $destriped 
    rm acolite_run*txt ${base}-acolite.log acolite.parameters polymer.parameters
    rm libenvironment-variables.so
fi

if $withtoolbox; then
    echo Progress[%%]: 100.0 : done
    echo $(pwd)/$newname
else
    echo $newname
fi
echo "done"
