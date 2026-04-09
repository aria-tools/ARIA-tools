# 1. Standard library imports
import warnings
from pathlib import Path
from typing import Optional, Tuple, Union

# 2. Third-party imports
import numpy as np
import rioxarray
import scipy.ndimage
from numpy.typing import NDArray
from osgeo import gdal, osr, gdal_array

#  READ/WRITE GDAL UTILITIES


def get_GUNW_attr(filename: Union[str, Path],
                  proj: Optional[str] = None,
                  xres: Optional[float] = None,
                  yres: Optional[float] = None,) -> dict:
    """
    Use GDAL to get raster metadata

    Parameters
    ----------
    filename : str
        path to raster

    Returns
    -------
    raster_attr : dict
        raster attribute dict
        [path, nodata, len, wid, snwe, lon_spacing, lat_spacing, projection]
    """

    warp_kwargs = dict(
        format="MEM",
        dstSRS=proj,
        multithread=False,
    )

    # set pixel spacing
    if xres is not None and yres is not None:
        warp_kwargs["xRes"] = xres
        warp_kwargs["yRes"] = yres

    # Use GDAL to read GUNW netcdf
    if proj is not None:
        ds = gdal.Warp("", str(filename), **warp_kwargs)
    else:
        ds = gdal.Open(filename, gdal.GA_ReadOnly)

    # Get GUNW Raster attributes
    nodata = ds.GetRasterBand(1).GetNoDataValue()

    # Get raster geographical information
    transform = ds.GetGeoTransform()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    snwe = [transform[3] + ysize * transform[5], transform[3],
            transform[0], transform[0] + xsize * transform[1]]
    lon_spacing = transform[1]
    lat_spacing = transform[5]

    projection = ds.GetProjection()

    # wrap raster info in dict
    raster_attr = {'PATH': filename,
                   'NODATA': nodata,
                   'LENGTH': ysize,
                   'WIDTH': xsize,
                   'SNWE': snwe,
                   'LON_SPACING': lon_spacing,
                   'LAT_SPACING': lat_spacing,
                   'PROJECTION': projection}

    # close
    ds = None

    return raster_attr


def get_GUNW_array(filename: Union[str, Path],
                   proj: str = "EPSG:4326",
                   nodata: Optional[float] = None,
                   subset: Optional[slice] = None,
                   xres: Optional[float] = None,
                   yres: Optional[float] = None,
                   align_to_grid: bool = False,
                   resample=gdal.GRA_NearestNeighbour,
                   as_xarray: bool = False,
                   varname: str = "connectedComponents",
                   mask: Optional[np.ndarray] = None,
                   ) -> np.ndarray:
    """
    Load a GUNW raster, optionally applying a binary mask, and reprojecting 
    to a consistent target grid. By default returns a NumPy array. If 
    `as_xarray=True`, returns an xarray.Dataset opened via the rasterio engine.
    """

    # Discover source nodata (if present)
    src = gdal.Open(str(filename), gdal.GA_ReadOnly)
    band = src.GetRasterBand(1)
    src_nodata = band.GetNoDataValue()

    # --- NEW: Apply the optional mask in native resolution before warping ---
    if mask is not None:
        # Create an in-memory copy of the source raster
        driver = gdal.GetDriverByName('MEM')
        warp_input = driver.CreateCopy('', src)
        masked_band = warp_input.GetRasterBand(1)
        arr = masked_band.ReadAsArray()

        # Validate dimensions
        if mask.shape != arr.shape:
            raise ValueError(f"Mask shape {mask.shape} does not match source shape {arr.shape}.")

        # Determine the safest NoData value to fill masked pixels with
        fill_val = src_nodata if src_nodata is not None else (nodata if nodata is not None else 0.0)

        # Apply mask: where mask == 1 (Valid), keep data; else replace with fill_val
        arr = np.where(mask == 1, arr, fill_val)

        # Write the cleanly masked array back into our MEM dataset
        masked_band.WriteArray(arr)
        
        # Ensure GDAL knows about our NoData value for the Warp step
        if src_nodata is None:
            masked_band.SetNoDataValue(fill_val)
            src_nodata = fill_val  # Update so warp_kwargs handles it correctly below
    else:
        # Standard operation: just point Warp to the file path
        warp_input = str(filename)

    src = None  # Safely release the original file lock

    # --- Original Warp and Processing Logic ---
    warp_kwargs = dict(
        format="MEM",
        dstSRS=proj,
        resampleAlg=resample,
        multithread=False,
        outputType=gdal.GDT_Float32
    )

    # set pixel spacing
    if xres is not None and yres is not None:
        warp_kwargs["xRes"] = xres
        warp_kwargs["yRes"] = yres
        if align_to_grid:
            warp_kwargs["targetAlignedPixels"] = True

    # Use explicit nodata handling for float rasters
    if nodata is not None:
        warp_kwargs["dstNodata"] = nodata
        warp_kwargs["srcNodata"] = src_nodata if src_nodata is not None else nodata

    # Reproject to target grid in-memory using our warp_input (either file path or MEM dataset)
    ds = gdal.Warp("", warp_input, **warp_kwargs)

    if not as_xarray:
        data = ds.ReadAsArray()
        ds = None
        if subset:
            data = data[subset]
        return np.ma.masked_equal(data, nodata) if nodata is not None else data

    # ---- xarray / rasterio path ----
    # Write to an in-memory GeoTIFF so rasterio can open it
    vsipath = "/vsimem/_gunw_tmp.tif"
    gdal.Translate(
        vsipath, ds, format="GTiff",
        creationOptions=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=IF_SAFER"]
    )
    ds = None  # release MEM dataset

    da = rioxarray.open_rasterio(vsipath, masked=True)

    # Optional subset (y, x) indexing; keep it simple if provided as a slice/tuple
    if subset:
        da = da.isel(y=subset[0], x=subset[1]) if isinstance(subset, tuple) else da[subset]

    # Squeeze single-band and name
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)

    # make sure it has a variable name
    da.name = varname

    # Ensure nodata encoded (use NaN if provided)
    if nodata is not None:
        da = da.rio.write_nodata(nodata)

    # Clean up the vsimem file
    gdal.Unlink(vsipath)

    return da.to_dataset(name=varname)


def write_GUNW_array(output_filename: Union[str, Path],
                     array: np.ndarray,
                     snwe: list,
                     nodata: Optional[str] = 'NAN',
                     format: Optional[str] = 'ENVI',
                     epsg: Optional[str] = 'EPSG:4326',
                     add_vrt: Optional[bool] = True,
                     verbose: Optional[bool] = False,
                     update_mode: Optional[bool] = True) -> None:
    """
    Use GDAL to write raster

    Parameters
    ----------
    output_filename : str
        path to raster
    array : array
        numpy array of raster to be written
    snwe  : list
        (South, North, West, East) bounds of raster to be written
    nodata : str
        value or nan for NODATA (used for VRT creation)
    format : str
        output raster format, default is ENVI
    epsg : str
        projection epsg, default is 'EPSG:4326' for WGS84
    add_vrt : bool
        flag to create VRT for output raster [True/False]
    verbose : bool
        print info messages [True/False]
    update_mode : bool
        flag to overwrite the existing files [True/False]
    """

    array_type = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)

    # Output path
    output = Path(output_filename).absolute()
    output_vrt = output.with_suffix('.vrt')

    if update_mode:
        [print(f'Remove {output}') if verbose else None]
        output.unlink(missing_ok=True)
        if add_vrt:
            [print(f'Remove {output_vrt}') if verbose else None]
            output_vrt.unlink(missing_ok=True)

    # Get lat, lon pixel spacing
    num_bands = 1
    if len(array.shape) > 2:
        num_bands = array.shape[0]
        x_step = (snwe[3] - snwe[2]) / array.shape[2]
        y_step = (snwe[0] - snwe[1]) / array.shape[1]
    else:
        x_step = (snwe[3] - snwe[2]) / array.shape[1]
        y_step = (snwe[0] - snwe[1]) / array.shape[0]

    # Geotransform
    geo = (snwe[2], x_step, 0, snwe[1], 0, y_step)
    srs = osr.SpatialReference()
    epsg_int = int(epsg.split(":")[-1])
    srs.ImportFromEPSG(epsg_int)  # set projection

    # Write
    driver = gdal.GetDriverByName(format)
    if len(array.shape) > 2:
        out_ds = driver.Create(
            str(output),
            array.shape[2],
            array.shape[1],
            array.shape[0],
            array_type)
    else:
        out_ds = driver.Create(
            str(output),
            array.shape[1],
            array.shape[0],
            1,
            array_type)

    out_ds.SetProjection(srs.ExportToWkt())
    out_ds.SetGeoTransform(geo)

    if verbose:
        print(f'Writing {output}')

    for i in range(num_bands):
        band = out_ds.GetRasterBand(i + 1)
        if num_bands > 1:
            band.WriteArray(array[i])
        else:
            band.WriteArray(array)
        band.FlushCache()
        band.ComputeStatistics(False)

    # Close
    out_ds = None

    if add_vrt:
        # Build virtual VRT
        vrt = gdal.BuildVRT(str(output_vrt), str(output), srcNodata=nodata)
        vrt.FlushCache()
        vrt = None


def snwe_to_extent(snwe: list) -> list:
    '''
    Convert SNWE to extent for matplotlib plotting
    '''
    extent = [snwe[2], snwe[3], snwe[0], snwe[1]]

    return extent


def _nan_filled_array(masked_array):
    masked_array.fill_value = np.nan
    return masked_array.filled()


def lalo2xy(lat: np.float32,
            lon: np.float32,
            data_snwe: list,
            latlon_step: list,
            rounding_method: Optional[str] = 'floor') \
        -> Tuple[gdal.GDT_Float32, gdal.GDT_Float32]:
    """
    Georeferenced coordinates to image space coordinates.
    GDAL raster starting point is the upper left corner.

    Parameters
    ----------
    lat : float
        search latitude
    lon : float
        search longitude
    data_snwe : list
        [South, North, West, East] bounds of raster
        North (y0) and West(x0) as reference point
    latlon_step : list
        pixel spacing [latitude_spacing, longitude_spacing]
    rounding_method : str
        rounding method, default is 'floor' other option is 'around'
        Read notes below. TODO. test different rounding routines

    Returns
    -------
    x : int
        coordinate in image space (x-axis/columns, direction of width)
        from x0 (top up column)
    y : int
        coordinate in image space (y-axis/rows, direction of length)
        from y0 (top left row)
    """

    # np.floor works better with points and raster - Need to check why
    # but with two rasters sometimes one pixel is missing or is redundant
    if rounding_method == 'floor':
        x = int(np.floor((lon - data_snwe[2]) / latlon_step[1] + 0.01))
        y = int(np.floor((lat - data_snwe[1]) / latlon_step[0] + 0.01))

    # np.around works better with two rasters
    # test it out, I think it has something to how numpy floor is
    # rounding negative values
    # example np.around(-125.2) = -125 np.floor(-125.2) = -126
    # np.around(125.6) = 126, np.floor(125.6) = 125
    elif rounding_method == 'around':
        x = int(np.around((lon - data_snwe[2]) / latlon_step[1] + 0.01))
        y = int(np.around((lat - data_snwe[1]) / latlon_step[0] + 0.01))

    return x, y


# Extract overlap bounds
def frame_overlap(snwe1: list,
                  snwe2: list,
                  latlon_step1: list,
                  latlon_step2: list,
                  latlon_step: Optional[list] = [-0.000833334, 0.000833334]
                  ) -> Tuple[tuple, tuple]:
    """
    Parameters
    ----------
    Use raster metadata to find overlap between two Images
    snwe1 : list
        [South, North, West, East] bounds of Image-1
    snwe2 : list
        [South, North, West, East] bounds of Image-2
    latlon_step1 : list
        latitude and longitude pixel spacing of Image-1
    latlon_step2 : list
        latitude and longitude pixel spacing of Image-2

    Returns
    -------
    subset1 : np.slice
        Intersection subset [y1:y2, x1:x2] for the Image-1
    subset2 : np.slice
        Intersection subset [y1:y2, x1:x2] for the Image-2
    """

    snwe = np.vstack([snwe1, snwe2])
    # Find overlap bounds
    overlap_snwe = [np.max(snwe[:, 0]), np.min(snwe[:, 1]),
                    np.max(snwe[:, 2]), np.min(snwe[:, 3])]

    # Georeferenced space to image coordinate space
    # Frame-1
    x1, y1 = lalo2xy(overlap_snwe[1], overlap_snwe[2], snwe1, latlon_step1)
    # Frame-2
    x2, y2 = lalo2xy(overlap_snwe[1], overlap_snwe[2], snwe2, latlon_step2)

    # Overlap bounds - force overlaps to have same dimensions
    # latlon_spacing sometimes diff at 13th decimal
    length = int(round((overlap_snwe[0] - overlap_snwe[1]) / latlon_step[0]))
    width = int(round((overlap_snwe[3] - overlap_snwe[2]) / latlon_step[1]))

    subset1 = np.s_[y1:y1 + length, x1:x1 + width]
    subset2 = np.s_[y2:y2 + length, x2:x2 + width]

    return subset1, subset2


def combine_data_to_single(data_list: list,
                           snwe_list: list,
                           latlon_step_list: list,
                           method: Optional[str] = 'mean',
                           latlon_step: Optional[list] =
                           [-0.000833334, 0.000833334]) \
        -> Tuple[NDArray, NDArray, list]:
    """
    Merge multiple arrays to one array. Combine them in ndarray, then apply
    function along the n_layers axis

    Parameters
    ----------
    data_list : list
        list of arrays containing raster values
    snwe_list : list
        list of arrays containing snwe (extent) values
    latlon_step_list : list
        list of arrays containing pixel spacing in lat and lon direction
        for each dataset
    method : str
        method to merge overlapping pixes, use mean, min, max etc..
        TODO: need to refine this part of code

    Returns
    -------
    comb_data : ndarray
        combined data [n_frames, length, width]
    SNWE : array
        extent of the combined data
    latlon_step : array
        pixel spacing in lat, lon of combined data

    """
    # Get the maximum extent of all data
    n = len(data_list)
    snwe_all = np.squeeze([snwe for snwe in snwe_list])

    SNWE = np.array([np.min(snwe_all[:, 0]), np.max(snwe_all[:, 1]),
                     np.min(snwe_all[:, 2]), np.max(snwe_all[:, 3])]).T

    length = abs(int(np.around((SNWE[1] - SNWE[0]) / latlon_step[0] + 0.01)))
    width = abs(int(np.around((SNWE[2] - SNWE[3]) / latlon_step[1] + 0.01)))

    # create combined data array
    # handle if 3D metadata layer
    if len(data_list[0].shape) > 2:
        comb_data = np.empty((n, data_list[0].shape[0], length, width),
                             dtype=np.float64) * np.nan
    else:
        comb_data = np.empty((n, length, width), dtype=np.float64) * np.nan
    for i, data in enumerate(data_list):
            x, y = np.abs(lalo2xy(SNWE[1], SNWE[2], snwe_list[i],
                                  latlon_step_list[i], 'around'))
            x, y = int(x), int(y)
            
            # handle if 3D metadata layer
            if len(data.shape) > 2:
                y_end = min(y + data.shape[1], comb_data.shape[2])
                x_end = min(x + data.shape[2], comb_data.shape[3])
                comb_data[
                    i, 0:data.shape[0], y:y_end, x:x_end
                ] = data[:, :y_end - y, :x_end - x]
            else:
                y_end = min(y + data.shape[0], comb_data.shape[1])
                x_end = min(x + data.shape[1], comb_data.shape[2])
                comb_data[
                    i, y:y_end, x:x_end
                ] = data[:y_end - y, :x_end - x]

    # Apply warning filters globally to the thread pool
    # instead of using a context manager
    # because Dask threads leak Python context managers and
    # cause the warning to bleed through.
    warnings.filterwarnings("ignore", message="Mean of empty slice")
    warnings.filterwarnings("ignore", message="All-NaN slice encountered")

    # combine using numpy
    if method == 'mean':
        comb_data = np.nanmean(comb_data, axis=0)
    elif method == 'median':
        comb_data = np.nanmedian(comb_data, axis=0)
    elif method == 'min':
        comb_data = np.nanmin(comb_data, axis=0)
    elif method == 'max':
        comb_data = np.nanmax(comb_data, axis=0)

    return comb_data, SNWE, latlon_step


def get_binary_nisar_mask_from_path(gdal_path):
    """
    Given a GDAL-style path, reconstructs the internal mask path.
    Uses spatial morphology to screen wide edge artifacts.
    """
    prefix, subpath = gdal_path.split('":')
    base, sep, rest = subpath.partition("unwrappedInterferogram")
    nisar_mask_path = f'{prefix}":{base}{sep}/mask'

    binary_mask = create_binary_nisar_mask(nisar_mask_path)

    pol = rest.strip("/").split("/")[0] if rest else "HH"
    unw_path = f'{prefix}":{base}{sep}/{pol}/unwrappedPhase'

    ds_unw = gdal.Open(unw_path, gdal.GA_ReadOnly)
    if ds_unw is not None:
        unw_arr = ds_unw.ReadAsArray()
        ds_unw = None

        # 1. Find the "core" of the artifact taper
        near_zeros = (np.abs(unw_arr) < 5e-5) & (unw_arr != 0)

        if np.any(near_zeros):
            # 2. Define a safe boundary containment zone (~75px)
            edge_zone = scipy.ndimage.binary_dilation(
                binary_mask == 0, iterations=75
            )
            core_artifacts = near_zeros & edge_zone

            # 3. Dilate the core to swallow the fading rest of the taper
            # which naturally exceeds the 5e-5 threshold
            full_artifacts = scipy.ndimage.binary_dilation(
                core_artifacts, iterations=40
            )

            # 4. Contain it within the edge zone and mask it
            full_artifacts = full_artifacts & edge_zone
            binary_mask = np.where(full_artifacts, 0, binary_mask)

    return binary_mask


def create_binary_nisar_mask(mask_path: str) -> np.ndarray:
    """
    Reads a 3-digit NISAR SAR mask file and decodes it into a
    binary mask in memory.
    
    1 = Valid (Secondary RSLC has data, ignores water status)
    0 = Invalid (Missing data in the secondary image ONLY, i.e., XX0)

    Args:
        mask_path (str): File path to the 3-digit mask file.

    Returns:
        np.ndarray: The decoded binary mask array.

    Raises:
        FileNotFoundError: If the mask file cannot be opened.
    """
    ds_mask = gdal.Open(mask_path, gdal.GA_ReadOnly)
    if not ds_mask:
        raise FileNotFoundError(f"Could not open mask file: {mask_path}")

    arr_mask = ds_mask.ReadAsArray()
    ds_mask = None  # Free GDAL dataset from memory

    # We only need to decode the least significant digit (Secondary RSLC)
    # Ignore the internal water mask for now
    sec_subswath = arr_mask % 10

    # Create and return binary mask
    # 1 if the last digit is not 0, otherwise 0
    binary_mask = np.where(sec_subswath != 0, 1, 0)

    return binary_mask


def apply_mask_and_write(
    vrt_unw_path: str,
    vrt_conn_path: str,
    binary_mask: np.ndarray,
    out_unw_path: str,
    out_conn_path: str,
    multiply_unw_by: int = 1
    ) -> Tuple[str, str]:
    """
    Applies a NISAR binary mask to VRT datasets, writes the results to temp
    files and OVERWRITES the original VRT files to point to the new temp files

    Args:
        vrt_unw_path (str): Path to the source unwrapped phase VRT to overwrite.
        vrt_conn_path (str): Path to the source conncomp VRT to overwrite.
        binary_mask (np.ndarray): The decoded binary mask array in memory.
        out_unw_path (str): Filepath to save the intermediate masked unw TIF.
        out_conn_path (str): Filepath to save the intermediate masked conn TIF.

    Returns:
        Tuple[str, str]: Paths to the newly overwritten VRT files.

    Raises:
        FileNotFoundError: If any input VRTs cannot be opened.
        ValueError: If the mask shape does not match the VRT shape.
    """
    ds_unw = gdal.Open(vrt_unw_path, gdal.GA_ReadOnly)
    ds_conn = gdal.Open(vrt_conn_path, gdal.GA_ReadOnly)

    if not all([ds_unw, ds_conn]):
        raise FileNotFoundError("One or more input VRTs could not be opened.")

    # Extract spatial metadata from the reference dataset
    geo_transform = ds_unw.GetGeoTransform()
    projection = ds_unw.GetProjection()
    cols = ds_unw.RasterXSize
    rows = ds_unw.RasterYSize

    if not (binary_mask.shape == (rows, cols)):
        raise ValueError("Binary mask shape does not match VRT dimensions.")

    # Apply the mask directly to the arrays read from the VRTs
    masked_unw = ds_unw.ReadAsArray() * binary_mask
    if multiply_unw_by == -1:
        masked_unw = masked_unw * -1
    masked_conn = ds_conn.ReadAsArray() * binary_mask

    # Driver for writing standard GeoTIFFs
    driver = gdal.GetDriverByName("GTiff")

    def _write_geotiff(
        out_path: str, data_array: np.ndarray, gdal_type: int
    ) -> None:
        """Helper to physically write the array to disk and safely close it."""
        out_ds = driver.Create(str(out_path), cols, rows, 1, gdal_type)
        out_ds.SetGeoTransform(geo_transform)
        out_ds.SetProjection(projection)
        
        band = out_ds.GetRasterBand(1)
        band.WriteArray(data_array)
        band.SetNoDataValue(0)
        band.FlushCache()
        
        out_ds = None  # Safely close the file lock

    # Write out the intermediate .tif files safely
    _write_geotiff(out_unw_path, masked_unw, ds_unw.GetRasterBand(1).DataType)
    _write_geotiff(out_conn_path, masked_conn, ds_conn.GetRasterBand(1).DataType)

    # CRITICAL: Free memory and release file locks on the original VRTs 
    # BEFORE we attempt to overwrite them in the next step.
    ds_unw = None
    ds_conn = None

    # Overwrite the original VRTs to point to our newly created TIFs
    gdal.BuildVRT(str(vrt_unw_path), str(out_unw_path))
    gdal.BuildVRT(str(vrt_conn_path), str(out_conn_path))

    return vrt_unw_path, vrt_conn_path
