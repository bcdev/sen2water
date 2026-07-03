# Sen2Water processor decomposition analysis

Sen2Water is composed of data processors that each have an internal structure. Purpose of this analysis is to identify intermediate data items, usually raster data with one or several bands, and respective units of functions that generate the intermediates from inputs or other intermediates. These intermediates are candidates for "breakpoints", data items that can be written into a file during processing. The units of functions are candidates for processing units of the EOPF framework or of lower level processing steps.

The scope of the analysis comprises

- msiresampling
- TGC
- Idepix
- C2RCC
- ACOLITE
- POLYMER
- s2wswitching

## Idepix

Idepix reads the L1C resampled to 60m and accesses respective tiles of a DEM and a static water mask and generates a flag band pixel_classif_flags.
Idepix for MSI is a SNAP operator that is composed of serveral internal operators. They are candidates for units. Their outputs are candidates for intermediates.

Candidate names for units are

| Unit (intermediate name) | Input bands                                               | Output bands                        | Computational service                | SNAP Operator                |
|--------------------------|-----------------------------------------------------------|-------------------------------------|--------------------------------------|------------------------------|
| l1c (product)            | -                                                         | B1,..,B12,sza,saa,vza,vaa,geocoding | main input                           | -                            | 
| dem (ADF)                | -                                                         | dem.elevation                       | ADF                                  | -                            | 
| watermask (ADF)          | -                                                         | watermask                           | ADF                                  | -                            |
| elevation                | geocoding, dem.elevation                                  | elevation                           | Interpolate DEM altitude to grid     | AddElevationOp               |
| slope_aspect_orientation | elevation, geocoding                                      | slope, aspect, orientation          | Determine slope, aspect, orientation | SlopeAspectOrientationOp     |
| idepix_classification    | B1,..,B12,sza,saa,vza,vaa,elevation,watermask             | raw.p_c_f                           | Classify by threshold tests          | S2ClassificationOp           | 
| idepix_filter_buffer     | raw.p_c_f,B2,B7,B8,B8A,B11                                | filtered.p_c_f                      | Filter spatially, apply buffer       | S2IdepixCloudPostProcessOp   | 
| idepix_mountain_shadow   | sza,saa,slope,aspect,orientation                          | mountain_shadow_flag                | Determine mountain shadow            | S2IdepixMountainShadowOp     | 
| idepix_cloud_statistics  | B8A,B3,geocoding,filtered.p_c_f,sza,saa,vza,vaa           | best_offset (triple)                | Determine cloud statistics           | S2IdepixPreCloudShadowOp     | 
| idepix_cloud_shadow      | sza,saa,vza,vaa                                           | sza_mean,saa_mean (scalars)         | Determine cloud shadow               | S2IdepixCloudShadowOp        | 
| idepix_shadow_clustering | B8A,B3,elevation,raw.p_c_f,best_offset,sza_mean,saa_mean  | shadow_flags (flag_band in SNAP Op) | Determine cloud shadow by clustering | S2IdepixPostCloudShadowOp    | 
| idepix_combine_flags     | filtered.p_c_f,mountain_shadow_flag,shadow_flag,raw.p_c_f | pixel_classif_flags                 | Combine flags                        | S2IdepixPostProcessOp        | 

There is an internal resampling for cloud shadow in case the input is not at 60m. In case the input is not in 60m there is a recursive call to MSI Idepix to get a classification at 60m. This case is not required for Sen2Water as long as it will only be used at 60m. The recursive call is not considered in this analysis.

## msiresampling

msiresmapling reads the L1C, resamples the reflectance bands to a common resolution (of 60m for Sen2Water), carefully aggregates angles from the band-specific and detector-specific angles of the L1C, and resamples ancillary data to the target grid.

The prototype implementation makes use of dask. It uses two classes Operator and Algorithm where Operator reads and creates xarray objects, Dataset for inputs and output, and dictionaries of DataArray for intermediates. The Algorithm is applied to dask arrays. Its implementation is a function on the level of numpy arrays applied to blocks. In msiresampling the same Algorithm is used in different places, and sometimes a step is implemented by alternative Algorithms depending on the source and target resolution. Several algorithms are applied per band separately. The top level operator creates the processing workflow.

- The `ResamplingOperator` is already on the level of an EOProcessingUnits. It can be rewritten as ResamplingUnit.
- There are two options for deconstruction:
- If we stay with Algorithms in this module we rewrite ResamplingOperator.run into ResamplingUnit.apply by changing the data type of parameters and return value to DataTree and the paths of variables within input and output datatree. It may also be necessary to transform input data dask arrays if they are read differently and not stacked.
- If we decompose msiresampling into many EOProcessingUnits then we wrap each algorithm with an EOProcessingUnit, so that the interface between units are xarray objects (xr.DataTrees) and unwrapping is handled within the unit. ResamplingUnit in this case orchestrates lower level units. We rewrite ResamplingOperator.run as well.
- The algorithms in MSIresampling can be used elsewhere. It is therefore preferable to keep them useful outside the EOPF framework and expose them as a library. The algorithms can then be wrapped as processing units and only the full ResamplingOperator's run function is completely reimplemented with the framework.

Candidate names for units are

| Unit (intermediate name)   | Input bands                                            | Output bands                                                      | Computational service                                                | Prototype class            |
|----------------------------|--------------------------------------------------------|-------------------------------------------------------------------|----------------------------------------------------------------------|----------------------------|
| MsiResampling              | B01,...,B12,B_detector_footprint_B1,...,               | B1,...,B12,B_detector_footprint_B1,...,                           | resampling to common target resolution                               | ResamplingOperator         |
|                            | quality_flags_B1,...,cloud_ice_flags,                  | quality_flags_B1,...,cloud_ice_flags,                             | considering detectors and angles                                     |                            |
|                            | msl,...,                                               | msl,...,msl_interpolated,...,                                     |                                                                      |                            |
|                            | sun_zenith,sun_azimuth,view_zenith,view_azimuth,       | sun_zenith,sun_azimuth,view_zenith_B1,...,view_azimuth_B1,...,    |                                                                      |                            |
|                            | B_ancillary_lost,...,spatial_ref_60m,...,y,x           | y,x,lat,lon,aux_latitude,aux_longitude,crs                        |                                                                      |                            |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| GetYX                      |                                                        | y,x                                                               | read data delayed (TBD whether this reading shall be delayed at all) | .values                    |
| GeoCoordinates             | y,x,crs                                                | lat,lon                                                           | transform UTM pixel coordinates to geo-coordinates                   | GeoCoordinates             |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| Downsampling               | B_detector_footprint_B1,...,master_detfoo              | resampled B_detector_footprint_B1,... or master_detfoo            | resample detector using majority, or find common detector,           | Downsampling               |
|                            |                                                        |                                                                   | per band, for higher resolutions                                     |                            |
| Upsampling                 | B_detector_footprint_B1,...                            | resampled B_detector_footprint_B1,...                             | resample detector, per band, for lower resolutions                   | Upsampling                 |
| MasterDetFoo               | master_detfoo,resampled B_detector_footprint_B1,...    | master_detfoo                                                     | mask out deviations in master_detfoo, for equal or lower resolutions | da.where                   |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| Downsampling               | B1,...B12,B_detector_footprint_B1,....,                | B1,...,B12                                                        | resample band using original and resampled detector for filtering,   | Downsampling               |
|                            | resampled B_detector_footprint_B1,... or master_detfoo |                                                                   | for higher resolutions                                               |                            |
| Upsampling                 | B1,...,B12                                             | B1,...,B12                                                        | resample band for lower resolutions                                  | Upsampling                 |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| Downsampling               | quality_flags_B1,...                                   | resampled quality_flags_B1,...                                    | resample band for higher resolution, using first, applied per band   | Downsampling               |
| Upsampling                 | quality_flags_B1,...                                   | resampled quality_flags_B1,...                                    | resample band for lower resolutions, using nearest, applied per band | Upsampling                 |
| Upsampling                 | cloud_ice_flags                                        | resampled cloud_ice_flags                                         | resample band for lower resolutions, using nearest, applied per band | Upsampling                 |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| GetAncillary               |                                                        | aux_latitude,aux_longitude,msl,...                                | read data delayed (TBD whether this reading shall be delayed at all) | .values                    |
| AncillaryInterpolation     | lat,lon,aux_latitude,aux_longitude,msl,...             | msl_interpolated,...                                              | interpolate from coarse geographic grid to target UTM grid           | AncillaryInterpolation     |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |
| TpInterpolation            | Bx,sun_zenith,sun_azimuth                              | sun_zenith, sun_azimuth                                           | interpolate angles to target grid                                    | TpInterpolation            | 
| GetViewingAngles           |                                                        | vza_B1,...,vaa_B1,...                                             | read data delayed (TBD whether this reading shall be delayed at all) | .values                    |
| ExpandViewingAngles        | vza_B1,...,vaa_B1,...                                  | extended vza_B1,...,vaa_B1,...                                    | extend angles per detector, delayed                                  | AnglesInterpolation.expand |
| AnglesInterpolation        | B_detector_footprint_B1,...,                           | view_zenith_B1,...,view_azimuth_B1,...                            | interpolate angles to target grid                                    | AnglesInterpolation        |
|                            | extended vza_B1,...,vaa_B1,...                         |                                                                   |                                                                      |                            |
| MeanAngles                 | view_zenith_B1,...,view_azimuth_B1,...                 | view_zenith_mean,view_azimuth_mean                                | mean angle for all bands                                             | MeanAngles                 |
| -------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------   |

## Polymer

- Processing functions are implemented using xarray. Therefore, `xarray.map_blocks` is used instead of `dask.array.map_blocks`
- We have to allow processing units with xarray based implementations if we want to use polymer functions without large rewrites (into numpy based functions)
- Or we analyse to which extend Polymer makes use of xarrray or to what extent it may even gain from using numpy.

## MSIDSF (reimplementation of parts of ACOLITE)

### Analysis

- launch_acolite.py is the entry point
- acolite_run.py collects settings and loops over inputs
- calls acolite_l1r format normalisation
  - identifies input type S2Resampling
  - reads lat, lon, vza, vaa, sza, saa
  - computes raa
  - reads reflectances, writes as rhot
  - reads ancillary from msl_interpolated etc. or msl etc. 
  - determines median (interpolated) or central (msl etc.)
- calls acolite_l2r for AC
  - determines proportion of invalid (blackfill) pixels from one (configured) band
  - optionally retrieves single aot value (optimise_aot_homogeneous)
  - optionally averages aot from spatio-temporal climatology
  - reads rsr lut, interpolates between wavelength
  - extracts meteo, interpolates, single value
  - reads and resamples DEM
  - calculates pressure from elevation
  - (selects par romix+rsky_t)
  - crop wind scalar to between 0.1 and 20
  - optionally crop sza and vza to limit suported by LUT
  - same as next with ozone lut
  - read Gas LUT dir/Gas/Gas_config.nc
  - interpolate transmittance value in LUT for pressure, wavelength, vza, sza
  - same with dir/LUB/WV/WV_config.nc
  - multiply transmittance
  - 
  - 