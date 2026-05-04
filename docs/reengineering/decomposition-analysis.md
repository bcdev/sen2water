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

Idepix for MSI is a SNAP operator that is composed of serveral internal operators. They are candidates for units. Their outputs are candidates for intermediates.

Candidate names for units are

| Unit (intermediate name) | Input bands                                               | Output Bands                        | Computational service                | SNAP Operator                |
|--------------------------|-----------------------------------------------------------|-------------------------------------|--------------------------------------|------------------------------|
| l1c (product)            | -                                                         | B1,..,B12,sza,saa,vza,vaa,geocoding | main input                           | -                            | 
| dem (ADF)                | -                                                         | dem.elevation                       | ADF                                  | -                            | 
| watermask (ADF)          | -                                                         | watermask                           | ADF                                  | -                            |
| elevation                | geocoding, dem.elevation                                  | elevation                           | Interpolate DEM altitude to grid     | AddElevationOp               |
| slope_aspect orientation | elevation, geocoding                                      | slope, aspect, orientation          | Determine slope, aspect, orientation | SlopeAspectOrientationOp     |
| idepix_classification    | B1,..,B12,sza,saa,vza,vaa,elevation,watermask             | raw.p_c_f                           | Classify by threshold tests          | S2ClassificationOp           | 
| idepix_filter_buffer     | raw.p_c_f,B2,B7,B8,B8A,B11                                | filtered.p_c_f                      | Filter spatially, apply buffer       | S2IdepixCloudPostProcessOp   | 
| idepix_mountain_shadow   | sza,saa,slope,aspect,orientation                          | mountain_shadow_flag                | Determine mountain shadow            | S2IdepixMountainShadowOp     | 
| idepix_cloud_statistics  | B8A,B3,geocoding,filtered.p_c_f,sza,saa,vza,vaa           | best_offset (triple)                | Determine cloud statistics           | S2IdepixPreCloudShadowOp     | 
| idepix_cloud_shadow      | sza,saa,vza,vaa                                           | sza_mean,saa_mean (scalars)         | Determine cloud shadow               | S2IdepixCloudShadowOp        | 
| idepix_shadow_clustering | B8A,B3,elevation,raw.p_c_f,best_offset,sza_mean,saa_mean  | shadow_flags (flag_band in SNAP Op) | Determine cloud shadow by clustering | S2IdepixPostCloudShadowOp    | 
| idepix_combine_flags     | filtered.p_c_f,mountain_shadow_flag,shadow_flag,raw.p_c_f | pixel_classif_flags                 | Combine flags                        | S2IdepixPostProcessOp        | 

There is an internal resampling for cloud shadow in case the input is not at 60m. In case the input is not in 60m there is a recursive call to MSI Idepix to get a classification at 60m. This case is not required for Sen2Water as long as it will only be used at 60m. The recursive call is not considered in this analysis.

