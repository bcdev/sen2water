## Version 0.2 (in development)

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
  (see SNAP issue ...)

## Version 0.1 (2024-01-09)

### Features

* Integration of ACOLITE, C2RCC, POLYMER
* Processing chain with S2Resampling, Idepix, all ACs, and merging
* Output product contains all AC results, yet un-merged
* Runtime structure with bin, etc, lib directory, sen2water.sh control script
* Separate directory with auxiliary data

