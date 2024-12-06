## Version 0.5.0 (2024-12-06)

### Fixes

### Enhancements

* Flag value OUT_OF_BOUNDS_SATURATED renamed to AC_OUT_OF_BOUNDS
* TOA Glint Correction added
* Switching made robust against missing bands in ACOLITE output due to low gas transmittance
* Windows support (with single-threaded POLYMER)
* Prepared for deployment as Sentinel Toolbox extension

## Version 0.4.2 (2024-09-24)

### Fixes

* Copernicus DEM handling for ante-meridian fixed in POLYMER (01KAB 20240103)
* Work-around for cfgrip reader of embedded ECMWF data with rotated longitude (01KAB), integrated into msiresampling and POLYMER, ticket https://github.com/ecmwf/cfgrib/issues/402
* Criterion for clouds flagged by C2RCC fixed, scale factor and offset applied (32MRU 20230808)
* s2wswitching main modified to avoid delayed NetCDF writing, attempt to fix sporadic infinite run.

### Enhancements

* MSI L1C reader updated to recognise the new namespace of N0511 in metadata XML files.

## Version 0.4.1 (2024-07-23)

### Fixes

* Copernicus DEM handling for certain latitudes (e.g. MSI granule 32VLM) patched, SNAP issue [SNAP-3752] raised
* libenvironment linked into working directory for Ubuntu 20

## Version 0.4 (2024-07-19)

### Enhancements

* POLYMER v4.7beta2 with configuration min_abs=-2
* ACOLITE 20240523 with configuration dsf_interface_reflectance=False
* msiresampling
* Global water mask, HR-OC mask for European coasts
* Switched back to conda because of GDAL dependency of ACOLITE

## Version 0.2 (continuous development)

### Enhancements

* Switching step added to select among AC outputs according 
  to water type and to static ocean-inland water mask
* S2Resampling step performance improved
* Conda Python environment exchanged with Python virtualenv 
  for smaller storage footprint
* Procedure defined how to build runtime package from source 
  and elements to be downloaded.

### Fixes

* S2Resampling exchanged against implementation that anchors  
  viewing angles grid at upper left corner of the image 
  (see SNAP issue ...) 
* Coordinates of CAMS and ECMWF data included in L1C added to 
  resampling intermediate to allow correct interpolation
  (see SNAP issue SNAP-3572)
* Switching step sometimes hangs. sequence in compute call changed.

## Version 0.1 (2024-01-09)

### Features

* Integration of ACOLITE, C2RCC, POLYMER
* Processing chain with S2Resampling, Idepix, all ACs, and merging
* Output product contains all AC results, yet un-merged
* Runtime structure with bin, etc, lib directory, sen2water.sh control script
* Separate directory with auxiliary data

