:: minimal gpt script for sen2water on Windows
:: adds idepix and c2rcc to classpath

@echo off
setlocal

set wd=%cd%
cd /d %~dp0%\..\..\..
set s2wdir=%cd%
cd /d %wd%

if "%blocksize%" == "" (
  set blocksize=610
  echo block size set to %blocksize%
) else (
  echo block size already set to %blocksize%
)

:: set PATH=%PATH%;%s2wdir%\lib\snap\snap\modules\lib\x86\jhdf5.dll
:: dir %s2wdir%\lib\snap\snap\modules\lib\x86\*

::echo "using java at %s2wdir%\lib\jre"
::echo "using snap at %s2wdir%\lib\snap"
::echo "working dir   %wd%"

::-XshowSettings:properties

%s2wdir%\lib\jre\bin\java ^
-cp %s2wdir%\lib\snap\snap\modules\*;%s2wdir%\lib\idepix\*;%s2wdir%\lib\c2rcc\* ^
-Duser.home=%wd% ^
-Djava.awt.headless=true ^
-Xmx6G ^
-Dsnap.mainClass=org.esa.snap.core.gpf.main.GPT ^
-Dsnap.dataio.reader.tileHeight=%blocksize% -Dsnap.dataio.reader.tileWidth=%blocksize% ^
-Dsnap.home=%s2wdir%\lib\snap ^
-Dsnap.userdir=%s2wdir% ^
-Dsnap.cachedir=%wd%\cache ^
-Dsnap.log.level=ERROR ^
org.esa.snap.runtime.Launcher %*
