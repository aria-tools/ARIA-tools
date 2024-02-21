import numpy as np
import warnings
from numpy.typing import NDArray

from typing import Optional, Tuple, Union
from osgeo import gdal, osr, gdal_array
from pathlib import Path

#  READ/WRITE GDAL UTILITIES


def get_GUNW_attr(filename: Union[str, Path]) -> dict:
    """
    Use GDAl to get raster metadata

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

    # Use GDAL to read GUNW netcdf
    ds = gdal.Open(filename)

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
                   nodata: Optional[float] = None,
                   subset: Optional[tuple] = None) -> NDArray:
    """
    Use GDAl to get raster data [as array]

    Parameters
    ----------
    filename : str
        path to raster
    subset : slice
        subset created with np.s_ TODO: insert tuple of (x1,x2, y1,y2)
        and then convert to slice with np.s_

    Returns
    -------
    data : array
        raster data 2D array [length x width]
    """

    # Use GDAL to read GUNW netcdf
    ds = gdal.Open(filename)

    data = ds.ReadAsArray()
    # close
    ds = None

    # Subset array
    if subset:
        data = data[subset]

    if nodata is not None:
        return np.ma.masked_equal(data, nodata)
    else:
        return data


def write_GUNW_array(output_filename: Union[str, Path],
                     array: np.ndarray,
                     snwe: list,
                     nodata: Optional[str] = 'NAN',
                     format: Optional[str] = 'ENVI',
                     epsg: Optional[int] = 4326,
                     add_vrt: Optional[bool] = True,
                     verbose: Optional[bool] = False,
                     update_mode: Optional[bool] = True) -> None:
    """
    Use GDAl to write raster

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
    epsg : int
        projection epsg, default is 4326 for WGS84
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
    srs.ImportFromEPSG(epsg)  # get projection

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
        -> Tuple[np.int16, np.int16]:
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
        # handle if 3D metadata layer
        if len(data.shape) > 2:
            comb_data[i, 0:data.shape[0], y:y + data.shape[1],
                      x: x + data.shape[2]] = data
        else:
            comb_data[i, y:y + data.shape[0], x: x + data.shape[1]] = data

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
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
