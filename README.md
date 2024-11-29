# Sen2Water

Sentinel-2 MSI pixel identification and atmospheric correction over ocean and inland water

## Input and output

Sen2Water processes Sentinel-2 MSI L1C granules with top-of-atmosphere reflectances. 
The input is expected in SAFE format, unpacked in a directory.
The 13 input bands of resolutions of 10m, 20m and 60m are processed into water-leaving 
reflectances over water areas at 60m resolution. In addition, flags are provided in 
three flag mask variables. The output data format is NetCDF4 with internal compression.

## Running the processor

The processor is called with

    . <path-to-mys2w>
    sen2water.sh <path-to-MSIL1C-SAFE-dir>

e.g.

    . /opt/sen2water-0.2/bin/mys2w
    sen2water.sh S2A_MSIL1C_20230611T103631_N0509_R008_T32UME_20230611T142059.SAFE

It will generate an output with the actual time stamp at the end of the file name:

    S2A_MSIL2W_20230611T103631_N0509_R008_T32UME_20240213T113100.nc

The runtime structure of the installation package is

    sen2water-0.5
    ├── bin
    │   ├── mys2w
    │   └── sen2water.sh
    │   └── sen2water.bat
    ├── etc
    ├── lib
    │   ├── acolite
    │   ├── c2rcc
    │   ├── conda
    │   ├── idepix
    │   ├── jre
    │   ├── msiresampling
    │   ├── polymer
    │   └── snap
    ├── auxdata
    ├── licenses
    └── README.md

## Installation

Sen2Water can be installed as stand-alone data processor. It is composed of elements that are 
downloaded from different sites and installed in an installation process. This makes sure that 
the individual licenses are recognised during installation. Ancillary data is provided as a 
separate package and installed either during installation or at runtime when needed.


