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

# parse parameters
input=
c2rccanc=
acoliteanc=
polymeranc=
dem=
withdetfoofilter=
chunksize=
withcleanup=false
while [ "$1" != "" ]; do
    if [ "$1" = "--c2rccanc" ]; then c2rccanc="$2"; shift 2
    elif [ "$1" = "--acoliteanc" ]; then acoliteanc="$2"; shift 2
    elif [ "$1" = "--polymeranc" ]; then polymeranc="$2"; shift 2
    elif [ "$1" = "--dem" ]; then dem="$2"; shift 2
    elif [ "$1" = "--withdetfoofilter" ]; then withdetfoofilter="$1"; shift 1
    elif [ "$1" = "--withcleanup" ]; then withcleanup=true; shift 1
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
    echo "--withdetfoofilter"
    echo "--withcleanup"
    echo "--chunksize 610 | 1830 | 915 | 366 | 305 | 183 | 122 | 61"
    exit 1
fi

# S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
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

echo "resampling to 60m ..."

time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py \
     --ancillary msl --ancillary tco3 --ancillary tcwv --ancillary u10 --ancillary v10 \
     $withdetfoofilter $chunksize $input $resampled

echo $resampled
echo
echo "Idepix cloud screening ..."

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
    $s2wdir/etc/idepix-graph.xml -Pdem="$dem" $resampled -t $idepix -f NetCDF4-BEAM

echo $idepix
echo
echo "C2RCC atmospheric correction ..."

if [ "$c2rccanc" == "embedded" ]; then
    useEcmwfAuxData=true
else
    useEcmwfAuxData=false
fi
time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
    -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
    $s2wdir/etc/c2rcc-graph.xml -Pdem="$dem" -PuseEcmwfAuxData="$useEcmwfAuxData" $resampled -t $c2rcc -f NetCDF4-BEAM

echo $c2rcc
echo
echo "ACOLITE atmospheric correction ..."

if [ "$acoliteanc" == "constant" ]; then
    s2auxiliarydefault=False
else
    s2auxiliarydefault=True
fi

cat $s2wdir/etc/acolite.parameters | sed \
    -e s,S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc,$resampled, \
    -e s,s2_auxiliary_default=True,s2_auxiliary_default=${s2auxiliarydefault}, \
    > acolite.parameters
time python -u $s2wdir/lib/acolite/launch_acolite.py --nogfx --cli --settings=acolite.parameters > ${base}-acolite.log

echo $acolite
echo
echo "POLYMER reformatting of cloud mask ..."

time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 4 -e \
    -Dsnap.dataio.reader.tileHeight=$blocksize -Dsnap.dataio.reader.tileWidth=$blocksize \
    $s2wdir/etc/polymer-mask-graph.xml $idepix -t $cloudmask -f NetCDF4-BEAM

echo $cloudmask
echo
echo "POLYMER atmospheric correction ..."

if [ "$polymeranc" == "" ]; then
    polymeranc=embedded
fi

cat $s2wdir/etc/polymer.parameters | sed -e s/S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_mask.nc/$cloudmask/ > polymer.parameters
rm -f $polymer
time python $s2wdir/lib/polymer/run-polymer.py polymer.parameters "$polymeranc" "$s2wdir/auxdata/dem/$dem" $resampled $polymer

echo $polymer
echo
echo "Sen2Water switching and output formatting ..."

time python -u $s2wdir/lib/msiresampling/sen2water/s2wswitching/main.py \
     $chunksize \
     $resampled $idepix $c2rcc $acolite $polymer $staticmask $s2w
newname=$(ncdump -h $s2w | grep ':id' | cut -d '"' -f 2)
mv $s2w $newname
s2w=$newname

if $withcleanup; then
  rm $resampled $idepix $c2rcc $acolite ${acolite/L2R/L1R} ${polymer/polymer/mask} $polymer
  rm acolite_run*txt ${resampled/resampled.nc/acolite.log} acolite.parameters polymer.parameters
fi

echo $s2w
echo "done"
