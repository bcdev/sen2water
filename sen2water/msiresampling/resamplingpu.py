from typing import Any, Optional, Literal, Tuple, Dict, Iterable, List

import numpy as np
import xarray as xr
import dask.array as da
from pyproj import Transformer, CRS

from eopf.computing.abstract import EOProcessingUnit
from eopf.computing.types import MappingAuxiliary, MappingDataType

from sen2water.msiresampling.ancillaryinterpolation import AncillaryInterpolation
from sen2water.msiresampling.anglesinterpolation import AnglesInterpolation
from sen2water.msiresampling.constants import MsiConstantsReengineering
from sen2water.msiresampling.downsampling import Downsampling
from sen2water.msiresampling.geocoordinates import GeoCoordinates
from sen2water.msiresampling.meanangles import MeanAngles
from sen2water.msiresampling.tpinterpolation import TpInterpolation
from sen2water.msiresampling.upsampling import Upsampling


class ResamplingPU(EOProcessingUnit):

    PROCESSOR_NAME = "MsiResampling"
    PROCESSOR_VERSION = "0.7.0"
    PROCESSOR_LEVEL = "Level-1 Product"
    PROCESSOR_MODEL = False

    INPUT_IDENTIFIER: str = "l1c"

    bands = MsiConstantsReengineering.bands
    resolutions = MsiConstantsReengineering.resolutions

    def __init__(self, identifier: Any = ""):
        super().__init__(identifier)
        self._update_history = True

    def run(
        self,
        inputs: MappingDataType,
        adfs: Optional[MappingAuxiliary] = None,
        mode: Optional[str] = None,
        *,
        resolution: int = 60,
        chunksize_in_meters = 36600,
        ancillary: Optional[List[str]]=None,
        downsampling: Literal['detectormean', 'first', 'min', 'max', 'mean', 'median'] = "detectormean",
        flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ] = "flagand",
        upsampling: Literal["nearest", "bilinear", "bicubic"] = "nearest",
        with_detfoo_filter: bool = False,
        merge_flags: bool = False,
        **kwargs
    ) -> MappingDataType:

        # check input

        if not "l1c" in inputs:
            raise KeyError("No 'l1c' in input of ResamplingMainPU.")
        l1c: xr.DataTree = inputs["l1c"]
        if not isinstance(l1c, xr.DataTree):
            raise TypeError("Input 'l1c' of ResamplingMainPU is not an xarray.DataTree.")

        # determine parameters

        overlap_depth = 2 if upsampling == 'bicubic' else 1 if upsampling == 'bilinear' else 0
        dims = {
            "y": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["y"],
            "x": l1c[f"measurements/reflectance/r{resolution}m"].coords.sizes["x"]
        }
        chunks = (chunksize_in_meters // resolution, chunksize_in_meters // resolution)

        # create result accumulator

        resampled = xr.DataTree(name="resampled")

        # add geo-coding with 1-d metric coordinates, CRS, 2-d geographic coordinates

        lat_data, lon_data = self._add_geocoding(l1c, resampled, ancillary, resolution, dims, chunks)

        # select detector per target pixel (with_detfoo_filter) or per target pixel per band

        if with_detfoo_filter:
            self._resample_detectors_with_detfoo_filter(l1c, resampled, resolution, chunksize_in_meters, dims, chunks)
        else:
            self._resample_detectors(l1c, resampled, resolution, chunksize_in_meters, dims, chunks)

        # resample reflectance bands B1 .. B12

        input_data_with_target_resolution = self._resample_reflectance(
            l1c, resampled,
            downsampling, upsampling, with_detfoo_filter, overlap_depth,
            resolution, chunksize_in_meters, dims, chunks,
        )

        # resample flag bands B_xxx_B1 .. B_zzz_B12

        self._resample_flags(
            l1c, resampled,
            flagdownsampling,
            resolution, chunksize_in_meters, dims, chunks,
        )

        self._resample_cloud_ice(
            l1c, resampled,
            flagdownsampling,
            resolution, chunksize_in_meters, dims, chunks,
        )

        # copy coarse grained meteo data

        resampled["conditions/meteorology/ecmwf"] = l1c["conditions/meteorology/ecmwf"].copy(inherit=False, deep=False)
        resampled["conditions/meteorology/cams"] = l1c["conditions/meteorology/cams"].copy(inherit=False, deep=False)

        # Removes mismatched encoding var when reading from zarr v2 and writing to zarr v3
        for node in resampled["conditions/meteorology"].subtree:
            for var in node.variables.values():
                var.encoding.clear()

        if ancillary:
            self._resample_ancillary(
                l1c, resampled,
                ancillary, lat_data, lon_data,
                dims, chunks,
            )
        del lat_data, lon_data

        self._resample_sun_angles(
            l1c, input_data_with_target_resolution, resampled,
            resolution, chunksize_in_meters, dims, chunks,
        )

        self._resample_viewing_angles(
            l1c, resampled,
            with_detfoo_filter,
            resolution, chunksize_in_meters, dims, chunks,
        )

        properties = l1c.attrs["stac_discovery"]["properties"]
        resampled.attrs = {
            "other_metadata": {
                "product_type": "S2_MSI_Level-1C",
                "platform": properties["platform"],
                "start_date": properties["start_datetime"],
                "stop_date": properties["end_datetime"],
            }
        }

        return {
            "resampled": resampled
        }

    def _add_geocoding(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            ancillary: list[str] | None,
            resolution: int,
            dims: dict[str, int],
            chunks: tuple[int, int],
    ) -> tuple[da.Array, da.Array]:

        epsg_code = l1c.attrs["stac_discovery"]["properties"]["proj:code"]
        # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
        # GT(1) w-e pixel resolution / pixel width.
        # GT(2) row rotation (typically zero).
        # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
        # GT(4) column rotation (typically zero).
        # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
        proj_transform = l1c.attrs["stac_discovery"]["properties"]["proj:transform"]

        geo_transform = " ".join(str(proj_transform[i]) for i in [2, 0, 1, 5, 3, 4])
        crs = CRS.from_authority(*epsg_code.split(":"))
        transformer = Transformer.from_crs(crs, CRS("EPSG:4326"))

        # The metric coordinates are numpy arrays, they are written as a single chunk.
        y_1d = l1c[f"measurements/reflectance/r{resolution}m/y"].data
        x_1d = l1c[f"measurements/reflectance/r{resolution}m/x"].data
        # Here we construct dask arrays to make the 2-d lat and lon chunked.
        y_data = da.from_array(y_1d, chunks=chunks[1])
        x_data = da.from_array(x_1d, chunks=chunks[0])
        y_grid = da.transpose(da.broadcast_to(y_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks[::-1]))
        x_grid = da.broadcast_to(x_data, (y_data.shape[0], x_data.shape[0]), chunks=chunks)

        lat_lon = GeoCoordinates().apply(
            x_grid,
            y_grid,
            transformer=transformer,
            new_axis=0,
            dtype=np.float64,
            chunks=(2, *chunks),  # TODO 2 is not optimal but the way the result knows that there are 2
        )
        lat_data = lat_lon[0]
        lon_data = lat_lon[1]

        y_da = xr.DataArray(
            y_1d,
            dims=["y"],
            attrs={
                "long_name": "y coordinate of projection",
                "standard_name": "projection_y_coordinate",
                "units": "m",
            },
        )
        x_da = xr.DataArray(
            x_1d,
            dims=["x"],
            attrs={
                "long_name": "x coordinate of projection",
                "standard_name": "projection_x_coordinate",
                "units": "m",
            },
        )
        crs_da = xr.DataArray(
            0, attrs={"crs_wkt": crs.to_wkt(), "GeoTransform": geo_transform}
        )
        lat_da = xr.DataArray(
            lat_data,
            dims=dims,
            attrs={"standard_name": "latitude",
                   "units": "degrees_north"},
        )
        lon_da = xr.DataArray(
            lon_data,
            dims=dims,
            attrs={"standard_name": "longitude",
                   "units": "degrees_east"},
        )

        resampled["measurements/reflectance/resampled/y"] = y_da
        resampled["measurements/reflectance/resampled/x"] = x_da
        resampled["measurements/reflectance/resampled/crs"] = crs_da
        resampled["measurements/reflectance/resampled/latitude"] = lat_da
        resampled["measurements/reflectance/resampled/longitude"] = lon_da
        resampled["conditions/mask/detector_footprint/resampled/y"] = y_da
        resampled["conditions/mask/detector_footprint/resampled/x"] = x_da
        resampled["conditions/mask/detector_footprint/resampled/crs"] = crs_da
        if ancillary:
            resampled["conditions/meteorology/resampled/y"] = y_da
            resampled["conditions/meteorology/resampled/x"] = x_da
            resampled["conditions/meteorology/resampled/crs"] = crs_da
        resampled["conditions/mask/l1c_classification/resampled/y"] = y_da
        resampled["conditions/mask/l1c_classification/resampled/x"] = x_da
        resampled["conditions/mask/l1c_classification/resampled/crs"] = crs_da
        resampled["conditions/geometry/resampled/y"] = y_da
        resampled["conditions/geometry/resampled/x"] = x_da
        resampled["conditions/geometry/resampled/crs"] = crs_da
        resampled["quality/mask/resampled/y"] = y_da
        resampled["quality/mask/resampled/x"] = x_da
        resampled["quality/mask/resampled/crs"] = crs_da

        return lat_data, lon_data


    def _resample_detectors(
        self,
        l1c: xr.DataTree,
        resampled: xr.DataTree,
        resolution: int,
        chunksize_in_meters: int,
        dims: Dict[str, int],
        chunks: tuple[int, int],
    ):
        for band in self.bands:
            band_resolution = self.resolutions[band]
            detector_footprint_band_name = (
                f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
            )
            band_chunksize = chunksize_in_meters // band_resolution
            detector_footprint_band_data = l1c[
                detector_footprint_band_name
            ].data.rechunk(band_chunksize)
            if resolution > band_resolution:
                factor = resolution // band_resolution
                resampled_detector = Downsampling().apply(
                    detector_footprint_band_data,
                    mode="majority",
                    factor=factor,
                    dtype=detector_footprint_band_data.dtype,
                    chunks=chunks,
                )

            elif resolution < band_resolution:
                factor = band_resolution // resolution
                resampled_detector = Upsampling().apply(
                    detector_footprint_band_data,
                    mode="nearest",
                    factor=(factor, factor),
                    src_image_shape=detector_footprint_band_data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=detector_footprint_band_data.dtype,
                    chunks=chunks,
                )
            else:
                resampled_detector = detector_footprint_band_data

            resampled[f"/conditions/mask/detector_footprint/resampled/{band}"] = (
                xr.DataArray(
                    resampled_detector,
                    dims=dims,
                    attrs={
                        "long_name": f"detector footprint of {band}",
                        "_FillValue": np.uint8(0),
                        "coordinates": "y x",
                        "grid_mapping": "crs",
                    },
                )
            )

    def _resample_detectors_with_detfoo_filter(
        self,
        l1c: xr.DataTree,
        resampled: xr.DataTree,
        resolution: int,
        chunksize_in_meters: int,
        dims: Dict[str, int],
        chunks: tuple[int, int],
    ):
        """Adds resampled detector index"""
        # The result shall be a single detfoo with fill value in the overlapping area without a common detfoo value.
        # Any of the source bands with the target resolution determines the master detector footprint value.
        # If any of the other bands does not have a contribution with this master detfoo value then the master detfoo pixel is set to invalid.
        master_band, master_detfoo = self._initialise_master_detfoo(l1c, resolution, chunks)

        for band in self.bands:
            if band == master_band:
                continue
            band_resolution = self.resolutions[band]
            detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
            band_chunksize = chunksize_in_meters // band_resolution
            detector_footprint_band_data = l1c[detector_footprint_band_name].data.rechunk(band_chunksize)
            if resolution > band_resolution:
                factor = resolution // band_resolution
                master_detfoo = Downsampling().apply(
                    detector_footprint_band_data,
                    master_detfoo,
                    mode="masterdetfoo",
                    factor=factor,
                    dtype=detector_footprint_band_data.dtype,
                    chunks=chunks
                )
            elif resolution < band_resolution:
                factor = band_resolution // resolution
                resampled_detector = Upsampling().apply(
                    detector_footprint_band_data,
                    mode="nearest",
                    factor=(factor, factor),
                    src_image_shape=detector_footprint_band_data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=detector_footprint_band_data.dtype,
                    chunks=chunks
                )
                master_detfoo = da.where(resampled_detector == master_detfoo, master_detfoo, 0)
            else:
                resampled_detector = detector_footprint_band_data
                master_detfoo = da.where(resampled_detector == master_detfoo, master_detfoo, 0)
        resampled["/conditions/mask/detector_footprint/resampled/master_detfoo"] = xr.DataArray(
            master_detfoo,
            dims=dims,
            attrs={
                "long_name": "master detector footprint",
                "_FillValue": np.uint8(0),
                "coordinates": "y x",
                "grid_mapping": "crs",
            }
        )

    def _initialise_master_detfoo(
            self, l1c: xr.DataTree, resolution: int, chunks: tuple[int, int]
    ) -> tuple[str, xr.DataTree]:
        master_band = None
        master_detfoo = None
        for band in self.bands:
            if self.resolutions[band] == resolution:
                master_band = band
                detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{resolution}m/{band}"
                master_detfoo_data = l1c[detector_footprint_band_name].data.rechunk(chunks)
                master_detfoo = xr.DataTree()
                master_detfoo['target_detfoo'] = xr.DataArray(
                    master_detfoo_data,
                    dims=l1c[detector_footprint_band_name].dims,
                    attrs={"long_name": "master detector footprint",
                           "_FillValue": np.uint8(0),
                           "coordinates": "y x",
                           "grid_mapping": "crs",
                    }
                )
                break
        if master_detfoo is None:
            raise KeyError(f"Missing input band with resolution {resolution} for master detfoo.")
        return master_band, master_detfoo

    def _resample_reflectance(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            downsampling: Literal[
                "detectormean", "mean", "median", "min", "max", "first",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
            ],
            upsampling: Literal["nearest", "bilinear", "bicubic"],
            with_detfoo_filter: bool,
            overlap_depth: int,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str,int],
            chunks: Tuple[int,int],
    ) -> da.Array:
        """Adds resampled reflectance bands, returns one input band with target resolution as dummy"""
        input_data_with_target_resolution: Optional[da.Array] = None
        for band in self.bands:
            band_resolution = self.resolutions[band]
            source_band = l1c[f"measurements/reflectance/r{band_resolution}m/{band}"]
            source_band_data = l1c[f"measurements/reflectance/r{band_resolution}m/{band}"].data.rechunk(chunksize_in_meters // band_resolution)
            factor = resolution // band_resolution
            if resolution > band_resolution:
                if downsampling == 'detectormean':
                    detector_footprint_band_name = f"/conditions/mask/detector_footprint/r{band_resolution}m/{band}"
                    detector_footprint_band_data = l1c[detector_footprint_band_name].data.rechunk(
                        chunksize_in_meters // band_resolution
                    )
                    resampled_band = Downsampling().apply(
                        source_band_data,
                        resampled["conditions/mask/detector_footprint/resampled/" + ("master_detfoo" if with_detfoo_filter else band)].data,
                        detector_footprint_band_data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=source_band.dtype,
                        chunks=chunks,
                    )
                else:
                    resampled_band = Downsampling().apply(
                        source_band_data,
                        mode=downsampling,
                        factor=factor,
                        is_reflectance=True,
                        dtype=source_band.dtype,
                        chunks=chunks,
                    )
            elif resolution < self.resolutions[band]:
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_band = Upsampling().apply(
                    source_band_data,
                    mode=upsampling,
                    factor=(factor, factor),
                    src_image_shape=source_band_data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    is_reflectance=True,
                    depth=overlap_depth,
                    trim=False,
                    dtype=source_band.dtype,
                    chunks=chunks
                )
            else:
                resampled_band = source_band_data
                if input_data_with_target_resolution is None:
                    input_data_with_target_resolution = source_band_data

            resampled[f"measurements/reflectance/resampled/{band}"] = xr.DataArray(
                resampled_band,
                dims=dims,
                attrs=source_band.attrs | {
                    "coordinates": "y x latitude longitude",
                    "grid_mapping": "crs",
                }
            )
        return input_data_with_target_resolution


    def _resample_ancillary(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            ancillary: Iterable[str],
            lat_data: da.Array,
            lon_data: da.Array,
            dims: Dict[str, int],
            chunks: Tuple[int,int],
    ):
        for anc_band_name in ancillary:
            anc_lat = l1c["conditions/meteorology/ecmwf/latitude"]
            anc_lon = l1c["conditions/meteorology/ecmwf/longitude"]
            anc_band = l1c[anc_band_name]
            anc_data = AncillaryInterpolation().apply(
                lat_data,
                lon_data,
                anc_lat=anc_lat.values,
                anc_lon=anc_lon.values,
                anc_data=anc_band.values,
                variable=anc_band,
                dtype=np.float32,
            )
            resampled[f"conditions/meteorology/resampled/{anc_band.name}"] = xr.DataArray(
                anc_data,
                dims=dims,
                attrs=anc_band.attrs.copy() | {
                    "coordinates": "y x",
                    "grid_mapping": "crs",
                }
            )

    def _resample_flags(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int,int],
    ):
        """Adds quality flags bands in target resolution"""
        for band in self.bands:
            band_resolution = self.resolutions[band]
            flag_band = f"quality/mask/r{band_resolution}m/{band}"
            flag_band_data = l1c[flag_band].data.rechunk(chunksize_in_meters // band_resolution)
            if resolution > band_resolution:
                factor = resolution // band_resolution
                resampled_flags = Downsampling().apply(
                    flag_band_data,
                    mode=flagdownsampling,
                    factor=factor,
                    dtype=flag_band_data.dtype,
                    chunks=chunks,
                )
            elif resolution < band_resolution:
                factor = band_resolution // resolution
                band_chunksize = chunksize_in_meters // band_resolution
                resampled_flags = Upsampling().apply(
                    flag_band_data,
                    mode='nearest',
                    factor=(factor, factor),
                    src_image_shape=flag_band_data.shape,
                    src_image_chunksize=(band_chunksize, band_chunksize),
                    depth=0,
                    dtype=flag_band_data.dtype,
                    chunks=chunks,
                )
            else:
                resampled_flags = flag_band_data

            resampled[f"quality/mask/resampled/{band}"] = xr.DataArray(
                resampled_flags,
                dims=dims,
                attrs=l1c[flag_band].attrs.copy() | {
                    "coordinates": "y x",
                    "grid_mapping": "crs",
                }
            )

    def _resample_cloud_ice(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            flagdownsampling: Literal[
                "mean", "median", "min", "max",
                "flagand", "flagor", "flagmedianand", "flagmedianor",
                "majority", "detectormean"
            ],
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int,int],
    ):
        """Adds cloud and ice flags in target resolution"""
        band = "conditions/mask/l1c_classification/r60m/b00"
        band_chunksize = chunksize_in_meters // 60
        flag_band_data = l1c[band].data.rechunk(band_chunksize)
        if resolution < 60:
            factor = 60 // resolution
            resampled_flags = Upsampling().apply(
                flag_band_data,
                mode='nearest',
                factor=(factor, factor),
                src_image_shape=flag_band_data.shape,
                src_image_chunksize=(band_chunksize, band_chunksize),
                depth=0,
                dtype=flag_band_data.dtype,
                chunks=chunks
            )
        else:
            resampled_flags = flag_band_data

        resampled[f"conditions/mask/l1c_classification/resampled/b00"] = xr.DataArray(
            resampled_flags,
            dims=dims,
            attrs=l1c[band].attrs.copy() | {
                    "coordinates": "y x",
                    "grid_mapping": "crs",
            }
        )

    def _resample_sun_angles(
            self,
            l1c: xr.DataTree,
            input_data_with_target_resolution: da.Array,
            resampled: xr.DataTree,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ):
        """Adds sun angles resampled from TP grid in target resolution"""

        angle_indices = {
            0: ("sza", "Sun Zenith Angle"),
            1: ("saa", "Sun Azimuth Angle"),
        }
        band = "conditions/geometry/sun_angles"
        for angle_idx, (angle_short_name, angle_long_name) in angle_indices.items():
            angle_data = l1c[band].data[angle_idx]
            band_chunksize = chunksize_in_meters // resolution
            resampled_angles = TpInterpolation().apply(
                input_data_with_target_resolution,
                tp_data=angle_data,
                resolution=resolution,
                tp_resolution=5000,
                image_shape=input_data_with_target_resolution.shape,
                image_chunksize=(band_chunksize, band_chunksize),
                dtype=angle_data.dtype,
            )
            resampled[f"conditions/geometry/resampled/{angle_short_name}"] = xr.DataArray(
                resampled_angles,
                dims=dims,
                attrs={
                    "long_name": angle_long_name,
                    "units": "degrees",
                    "_FillValue": np.nan,
                    "coordinates": "y x",
                    "grid_mapping": "crs",
                }
            )

    def _resample_viewing_angles(
            self,
            l1c: xr.DataTree,
            resampled: xr.DataTree,
            with_detfoo_filter: bool,
            resolution: int,
            chunksize_in_meters: int,
            dims: Dict[str, int],
            chunks: Tuple[int, int],
    ):
        """Resamples viewing angles per detector and adds viewing angles per band"""
        vza_accu = []
        vaa_accu = []

        view_angle_band_name = "conditions/geometry/viewing_incidence_angles"
        view_angle_band = l1c[view_angle_band_name]

        # Shape of the view_angle_band
        # "coordinates": [
        #   "detector",
        #   "angle",
        #   "y",
        #   "x",
        #   "band"
        # ],
        angle_indices = {
            "vza": 0,
            "vaa": 1
        }
        detectors = np.array(
            [int(d.removeprefix("d")) for d in view_angle_band.coords["detector"].values]
        )

        for band in self.bands:
            detector_footprint_band_name = "conditions/mask/detector_footprint/resampled/" + ("master_detfoo" if with_detfoo_filter else band)
            detector_footprint = resampled[detector_footprint_band_name].data
            extended_vza, extended_vaa = AnglesInterpolation().expand_angles_per_detector(
                view_angle_band.isel(angle=angle_indices["vza"]).sel(band=band).values,
                view_angle_band.isel(angle=angle_indices["vaa"]).sel(band=band).values,
            )
            resampled_vza = AnglesInterpolation().apply(
                detector_footprint,
                detectors=detectors,
                detector_angles=extended_vza,
                resolution=resolution,
                angles_resolution=5000,
                band=band,
                image_chunksize=chunks,
                dtype=np.float32,  # TBC
            )
            resampled_vaa = AnglesInterpolation().apply(
                detector_footprint,
                detectors=detectors,
                detector_angles=extended_vaa,
                resolution=resolution,
                angles_resolution=5000,
                band=band,
                is_azimuth_angle=True,
                image_chunksize=chunks,
                dtype=np.float32,  # TBC
            )
            # Shape of the view_angle_band
            # "coordinates": [
            #   "detector",
            #   "angle",
            #   "y",
            #   "x",
            #   "band"
            # ],
            resampled[f"conditions/geometry/resampled/vza_{band}"] = xr.DataArray(
                resampled_vza,
                dims=dims,
                attrs={
                    "long_name": f"Viewing incidence zenith angle for band {band}",
            "units": "degrees",
                    "_FillValue": np.nan,
                    "coordinates": "y x",
                    "grid_mapping": "crs",
                }
            )
            resampled[f"conditions/geometry/resampled/vaa_{band}"] = xr.DataArray(
                resampled_vaa,
                dims=dims,
                attrs={
                    "long_name": f"Viewing incidence azimuth angle for band {band}",
                    "units": "degrees",
                    "_FillValue": np.nan,
                    "coordinates": "y x",
                    "grid_mapping": "crs",
                }
            )
            vza_accu.append(resampled_vza)
            vaa_accu.append(resampled_vaa)

        vza_mean = MeanAngles().apply(*vza_accu, dtype=np.float32)
        vaa_mean = MeanAngles().apply(*vaa_accu, is_azimuth_angle=True, dtype=np.float32)

        resampled[f"conditions/geometry/resampled/vza_mean"] = xr.DataArray(
            vza_mean,
            dims=dims,
            attrs={
                "long_name": "Mean viewing incidence zenith angle",
                "units": "degrees",
                "_FillValue": np.nan,
                "coordinates": "y x",
                "grid_mapping": "crs",
            }
        )
        resampled[f"conditions/geometry/resampled/vaa_mean"] = xr.DataArray(
            vaa_mean,
            dims=dims,
            attrs={
                "long_name": "Mean viewing incidence azimuth angle",
                "units": "degrees",
                "_FillValue": np.nan,
                "coordinates": "y x",
                "grid_mapping": "crs",
            }
        )
