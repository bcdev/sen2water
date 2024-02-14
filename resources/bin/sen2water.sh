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
    echo "sen2water.sh <l1cpath>"
    echo "e.g."
    echo "sen2water.sh S2A_MSIL1C_20240104T103431_N0510_R108_T32UME_20240104T123149.SAFE"
    exit 1
fi

# S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
input=$1
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
staticmask=${s2wdir}/auxdata/s2w-mask/${granule:0:2}/s2w-mask-${granule}.tif
s2w=${base}-s2w.nc

echo "resampling to 60m ..."

#ln -sf $s2wdir/lib/snap/snap/modules/lib/amd64/libenvironment-variables.so
##    -Dsnap.production.concurrent=false -Dsnap.gpf.disableTileCache=true \
#time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -c 2048M -q 4 -e \
#    -Dsnap.gpf.disableTileCache=true \
#    $s2wdir/etc/s2resampling-graph.xml $input/MTD_MSIL1C.xml -t $resampled -f NetCDF4-BEAM
#rm libenvironment-variables.so

export PYTHONPATH=$s2wdir/lib/msiresampling:$PYTHONPATH
time python -u $s2wdir/lib/msiresampling/sen2water/msiresampling/main.py $input $resampled

echo $resampled
echo
echo "Idepix cloud screening ..."

time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 3 -e \
    $s2wdir/etc/idepix-graph.xml $resampled -t $idepix -f NetCDF4-BEAM

echo $idepix
echo
echo "C2RCC atmospheric correction ..."

time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 3 -e \
    $s2wdir/etc/c2rcc-graph.xml $resampled -t $c2rcc -f NetCDF4-BEAM

echo $c2rcc
echo
echo "ACOLITE atmospheric correction ..."

cat $s2wdir/etc/acolite.parameters | sed -e s/S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc/$resampled/ > acolite.parameters
time python -u $s2wdir/lib/acolite/launch_acolite.py --nogfx --cli --settings=acolite.parameters > ${base}-acolite.log

echo $acolite
echo
echo "POLYMER reformatting of cloud mask ..."

time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 3 -e \
    $s2wdir/etc/polymer-mask-graph.xml $idepix -t $cloudmask -f NetCDF4-BEAM

echo $cloudmask
echo
echo "POLYMER atmospheric correction ..."

cat $s2wdir/etc/polymer.parameters | sed -e s/S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_mask.nc/$cloudmask/ > polymer.parameters
#export PYTHONPATH=/home/yarn/appcache/sen2water2/sen2water-0.1/lib/polymer:$PYTHONPATH
rm -f $polymer
time python $s2wdir/lib/polymer/run-polymer.py polymer.parameters $resampled $polymer

echo $polymer
echo
echo "Sen2Water switching and output formatting ..."

if [ -e ${staticmask} ]; then
#    time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 3 -e \
#         $s2wdir/etc/s2w-merge-graph.xml \
#         -Sresampled=$resampled \
#         -Sidepix=$idepix \
#         -Sc2rcc=$c2rcc \
#         -Sacolite=$acolite \
#         -Spolymer=$polymer \
#         -Smask=${staticmask} \
#         -t $s2w -f NetCDF4-BEAM
    export PYTHONPATH=$s2wdir/lib/s2wswitching:$PYTHONPATH
    time python -u $s2wdir/lib/s2wswitching/sen2water/s2wswitching/main.py --copyinputs $resampled $idepix $c2rcc $acolite $polymer $staticmask $s2w
    #time python $s2wdir/lib/switching/main.py --copyinputs $resampled $idepix $c2rcc $acolite $polymer $staticmask $s2w
    newname=$(ncdump -h $s2w | grep ':id' | cut -d '"' -f 2)
    mv $s2w $newname
    s2w=$newname
else
    time gpt -J-Xmx6G -Dsnap.userdir=${s2wdir} -Dsnap.cachedir=$(pwd)/.snap/var -Dsnap.log.level=ERROR -c 4096M -q 3 -e \
         $s2wdir/etc/s2w-merge-graph.xml \
         -Sresampled=$resampled \
         -Sidepix=$idepix \
         -Sc2rcc=$c2rcc \
         -Sacolite=$acolite \
         -Spolymer=$polymer \
         -t $s2w -f NetCDF4-BEAM
fi

echo $s2w
echo "done"
