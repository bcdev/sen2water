:: sen2sater.bat Windows script to run the Sen2Water processing chain
:: that generates an MSI L2W from an MSI L1C.

@echo off
setlocal enableDelayedExpansion enableExtensions

set wd1=%cd%
cd %~dp0..%
set s2wdir=%cd%
cd %wd1%

echo Sen2Water 0.6.1 (%s2wdir%)

where /q sen2water.bat > Nul
if %ERRORLEVEL% neq 0 (
    echo setting up environment ...
    set PYTHONPATH=%s2wdir%\lib\msiresampling
    set PATH=%s2wdir%\bin;%s2wdir%\lib\snap\bin;%s2wdir%\lib\snap\snap\modules\lib\amd64;%PATH%
    set ECCODES_DEFINITION_PATH="%s2wdir%\lib\conda\share\eccodes\definitions"
    call %s2wdir%\lib\conda\condabin\activate.bat
)

:: parse parameters
set input=
set c2rccanc=
set acoliteanc=
set polymeranc=
set dem=
set withouttgc=
set withdetfoofilter=
set chunksize=
set withoutcleanup=
set withtoolbox=false
set outputdir=

:loop
if "%1" NEQ "" (
   if "%1" == "--c2rccanc" (
       set c2rccanc=%2
       shift
       shift
   ) else ( if "%1" == "--acoliteanc" (
       set acoliteanc=%2
       shift
       shift
   ) else ( if "%1" == "--polymeranc" (
       set polymeranc=%2
       shift
       shift
   ) else ( if "%1" == "--dem" (
       set dem=%2
       shift
       shift
   ) else ( if "%1" == "--withouttgc" (
       shift
       set withouttgc=true
   ) else ( if "%1" == "--withdetfoofilter" (
       shift
       set withdetfoofilter=true
   ) else ( if "%1" == "--withoutcleanup" (
       shift
       set withoutcleanup=true
   ) else ( if "%1" == "--withtoolbox" (
       shift
       set withtoolbox=true
   ) else ( if "%1" == "--outputdir" (
       set outputdir=%2
       shift
       shift
   ) else ( if "%1" == "--chunksize" (
       set chunksize=%2
       shift
       shift
   ) else ( if "%1" == "--help" (
       shift
   ) else ( if "%1:0:1%" == "-" (
       echo unknown parameter %1
       exit /b 1
   ) else ( if "!input!" == "" (
       set input=%~1
       shift
   ) else (
       echo unknown parameter %1
       exit /b 1
   )))))))))))))
   goto loop
)

if "!input!" == "" (
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
    exit /b 1
)

:: S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112.SAFE
if "!input:~-1!" == "/" (
    set input=!input:~0,-1!"
)
if "!input:~-4!" == ".xml" (
    call :dirname "!input!" input
)
call :basename "!input!" base0
if "%base0:~-5%" == ".SAFE" (
    set base=%base0:~0,-5%
) else (
    set base=%base0%
)
set dummy=dummytoavoiderrormessage

set granule=%base:~39,5%
set zone=%base:~39,2%
set p=%base:~0,3%
set y=%base:~11,4%
set m=%base:~15,2%
set d=%base:~17,2%
set H=%base:~20,2%
set I=%base:~22,2%
set S=%base:~24,2%
set resampled=%base%-resampled.nc
set destriped=%base%-resampled_TGC.nc
set idepix=%base%-idepix.nc
set c2rcc=%base%-c2rcc.nc
set acolite=%p%_MSI_%y%_%m%_%d%_%H%_%I%_%S%_S2R_L2R.nc
set cloudmask=%base%-mask.nc
set polymer=%base%-polymer.nc
set s2w=%base%-s2w.nc

set staticmask="%s2wdir%\auxdata\s2w-mask\%zone%\s2w-mask-%granule%.tif"
if not exist !staticmask! (
    set staticmask="%s2wdir%\auxdata\s2w-global-mask\%zone%\s2w-globalmask-%granule%.tif"
    if not exist !staticmask! (
        echo missing mask !staticmask!
        echo assuming ocean
        set staticmask=ocean
    )
)

if "!withtoolbox!" == "true" (
    if "!outputdir!" == "" (
        call :dirname !input! outputdir
    )
    cd !outputdir!
) else (
    set outputdir=%cd%
)
echo working directory !outputdir!

:: check whether we modify the software directory
echo.%cd% | find "%s2wdir%" > Nul && ( 
    echo "%~f0%" must not be started from within software dir "%s2wdir%"
    echo cd to some working directory outside the software installation, please
    exit /b 1
)

::ln -sf ${s2wdir}/lib/snap/snap/modules/lib/amd64/libenvironment-variables.so .

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 5.0 : resampling to 60m ...
) else (
    echo resampling to 60m ...
)

python -u %s2wdir%\lib\msiresampling\sen2water\msiresampling\main.py ^
     --ancillary msl --ancillary tco3 --ancillary tcwv --ancillary u10 --ancillary v10 ^
     !withdetfoofilter! !chunksize! !input! %resampled%

echo %resampled%
echo:

if "!withouttgc!" == "true" (
    set destriped=%resampled%
) else (
    if "!withtoolbox!" == "true" (
        echo Progress[%%]: 25.0 : TOA glint correction ...
    ) else (
        echo TOA glint correction ...
    )
    python %s2wdir%\lib\acolite\tgcparameters.py %input%
    python -u %s2wdir%\lib\acolite\tgc2.py --acolite_path %s2wdir%\lib\acolite --input %resampled% --output . ^
           --glint_threshold 0.0 --scaling True --aot_min 0.1 --process_high_sza False --estimate False ^
           --grid_write True --grid_files %base%-TGC-parameters.json --verbosity 2
    if exist %destriped% (
        echo %destriped%
        echo:
    ) else (
        echo "*** no %destriped% found, using resampled instead"
        set destriped=%resampled%
    )
)

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 45.0 : ACOLITE atmospheric correction ...
) else (
    echo ACOLITE atmospheric correction ...
)

if "!acoliteanc!" == "constant" (
    set s2auxiliarydefault=False
) else (
    set s2auxiliarydefault=True
)

powershell -Command "(gc %s2wdir%\etc\acolite.parameters) -replace 'S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc', '!destriped!' -replace 's2_auxiliary_default=True', 's2_auxiliary_default=!s2auxiliarydefault!' | Out-File -encoding ASCII acolite.parameters"

python -u %s2wdir%\lib\acolite\launch_acolite.py --nogfx --cli --settings=acolite.parameters > %base%-acolite.log

if not exist %acolite% (
    echo *** ACOLITE failed, retrying without TGC ...
    powershell -Command "(gc %s2wdir%\etc\acolite.parameters) -replace 'S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_60m.nc', '!resampled!' -replace 's2_auxiliary_default=True', 's2_auxiliary_default=!s2auxiliarydefault!' | Out-File -encoding ASCII acolite.parameters"
    python -u %s2wdir%\lib\acolite\launch_acolite.py --nogfx --cli --settings=acolite.parameters > %base%-acolite.log
    if not exist %acolite% (
        echo *** ACOLITE failed
        exit /b 1
    )
    set destriped=%resampled%
)

echo %acolite%
echo:

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 55.0 : Idepix cloud screening ...
) else (
    echo Idepix cloud screening ...
)

if "!dem!" == "" (
    set dem="Copernicus 90m Global DEM"
    echo dem is !dem!
)
if "!chunksize!" == "" (
    set blocksize=610
) else (
    for /f "tokens=2" %%G IN ('!chunksize!') DO set blocksize=%%G
)

:: -J-Xmx6G -Dsnap.userdir=%s2wdir% -Dsnap.cachedir=%cd%\.snap\var -Dsnap.log.level=ERROR

call gpt.bat -Dsnap.dataio.reader.tileHeight=%blocksize% -Dsnap.dataio.reader.tileWidth=%blocksize% ^
    -c 4096M -q 4 -e ^
    %s2wdir%\etc\idepix-graph.xml -Pdem=!dem! !destriped! ^
    -t %idepix% -f NetCDF4-BEAM

echo %idepix%
echo:

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 65.0 : C2RCC atmospheric correction ...
) else (
    echo C2RCC atmospheric correction ...
)

if "!c2rccanc!" == "embedded" (
    set useEcmwfAuxData=true
) else (
    set useEcmwfAuxData=false
)
call gpt.bat -Dsnap.dataio.reader.tileHeight=%blocksize% -Dsnap.dataio.reader.tileWidth=%blocksize% ^
    -c 4096M -q 4 -e ^
    %s2wdir%\etc\c2rcc-graph.xml -Pdem=!dem! -PuseEcmwfAuxData=!useEcmwfAuxData! !destriped! ^
    -t %c2rcc% -f NetCDF4-BEAM

echo %c2rcc%
echo:

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 75.0 : POLYMER atmospheric correction ...
) else (
    echo POLYMER atmospheric correction ...
)

call gpt.bat -Dsnap.dataio.reader.tileHeight=%blocksize% -Dsnap.dataio.reader.tileWidth=%blocksize% ^
   -c 4096M -q 4 -e ^
   %s2wdir%\etc\polymer-mask-graph.xml %idepix% ^
   -t %cloudmask% -f NetCDF4-BEAM

echo %cloudmask%

cd %s2wdir%\lib\polymer
python setup.py build_ext --inplace
cd !outputdir!

if "!polymeranc!" == "" (
    set polymeranc=embedded
)
powershell -Command "(gc %s2wdir%\etc\polymer.parameters.windows) -replace 'S2A_MSIL1C_20230929T103821_N0509_R008_T32UME_20230929T141112_mask.nc', '%cloudmask%' | Out-File -encoding ASCII polymer.parameters"

if exist %polymer% del %polymer%
python -u %s2wdir%\lib\polymer\run-polymer.py polymer.parameters !polymeranc! "%s2wdir%\auxdata\dem\!dem:~1,-1!" %resampled% %polymer%

echo %polymer%
echo:

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 85.0 : Sen2Water switching and output formatting ...
) else (
    echo Sen2Water switching and output formatting ...
)

python -u %s2wdir%\lib\msiresampling\sen2water\s2wswitching\main.py ^
     !chunksize! ^
     !destriped! %idepix% %c2rcc% %acolite% %polymer% !staticmask! %s2w%

if "!withtoolbox!" == "true" (
    echo Progress[%%]: 95.0 : renaming output ...
) else (
    echo renaming output ...
)

for /f "tokens=* USEBACKQ" %%f in (`ncdump -h %s2w% ^| find ":id"`) do (
    set "fid=%%f"
)
:: set "newname=%fid:~7,-3%" # netcdf4 instead of h5netcdf
set "newname=%fid:~14,-3%"
ren %s2w% %newname%

if "!withtoolbox!" == "true" (
    mkdir %temp%\s2w-output
    del %temp%\s2w-output\S2*
    copy %newname% %temp%\s2w-output
)

if not "!withoutcleanup!" == "true" (
    del %resampled% %idepix% %c2rcc% %acolite% %acolite:L2R=L1R% %polymer:polymer=mask% %polymer%
    del %base%-TGC-parameters.json !destriped! 
    del acolite_run*txt %base%-acolite.log acolite.parameters polymer.parameters
)

echo !newname!
if "!withtoolbox!" == "true" (
    echo Progress[%%]: 100.0 : done
    echo %cd%\%newname%
) else (
    echo %newname%
)
echo done

exit /b 0

:dirname p v
    setlocal enableextensions enabledelayedexpansion
    set _dir=%~dp1
    set _dir=%_dir:~0,-1%
    endlocal & set %2=%_dir%
    exit /b 0

:basename p v
    setlocal enableextensions enabledelayedexpansion
    set _dir=%~n1
    endlocal & set %2=%_dir%
    exit /b 0
