# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import logging
import math

import numpy as np
import scipy.interpolate

LOGGER = logging.getLogger(__name__)


def _compute_dem_range(dem_ds):
    """Compute the (min, max) elevation of a DEM dataset, skipping nodata.

    GDAL's ``ComputeRasterMinMax`` can silently include nodata pixels
    when the nodata value cannot be represented in the band's data type
    (e.g. NaN nodata on an Int16 band — former NaN pixels become 0 and
    are counted as valid).  This function reads the band array and masks
    nodata explicitly with numpy to ensure correct results.

    Parameters
    ----------
    dem_ds : gdal.Dataset
        Opened GDAL dataset for the DEM.

    Returns
    -------
    tuple of (float, float)
        ``(dem_min, dem_max)`` excluding nodata pixels.

    Raises
    ------
    ValueError
        If all DEM pixels are nodata.
    """
    band = dem_ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float64)
    nodata = band.GetNoDataValue()

    if nodata is not None and math.isnan(nodata):
        mask = ~np.isnan(arr)
    elif nodata is not None:
        mask = arr != nodata
    else:
        mask = np.ones(arr.shape, dtype=bool)

    if not np.any(mask):
        raise ValueError('All DEM pixels are nodata')

    dem_min = float(np.min(arr[mask]))
    dem_max = float(np.max(arr[mask]))
    LOGGER.debug('DEM range (nodata-aware): %.1f to %.1f m', dem_min, dem_max)
    return dem_min, dem_max


def _get_height_subset_indices(heights, dem_min, dem_max, pad=0):
    """Determine which height layers span the DEM elevation range.

    Parameters
    ----------
    heights : np.ndarray
        1-D array of height levels (ascending or descending).
    dem_min, dem_max : float
        Minimum / maximum DEM elevation.
    pad : int, optional
        Extra layers on each side of the bracket (default 0).

    Returns
    -------
    np.ndarray
        Integer indices into *heights* for the needed subset.
    """
    n = len(heights)
    if n <= 1:
        return np.arange(n)

    ascending = heights[-1] > heights[0]
    h_sorted = heights if ascending else heights[::-1]

    idx_lo = int(np.searchsorted(h_sorted, dem_min, side='right')) - 1
    idx_hi = int(np.searchsorted(h_sorted, dem_max, side='left'))

    idx_lo = max(0, idx_lo - pad)
    idx_hi = min(n - 1, idx_hi + pad)

    if not ascending:
        idx_lo, idx_hi = n - 1 - idx_hi, n - 1 - idx_lo

    return np.arange(idx_lo, idx_hi + 1)


class InterpCube(object):
    """Class to interpolate intersection of cube with DEM."""

    def __init__(self, inobj, hgtobj, latobj, lonobj, dem_range=None):
        """Init with h5py dataset.

        Parameters
        ----------
        inobj : array-like
            3-D data of shape (heights, lats, lons).
        hgtobj : array-like
            1-D height levels.
        latobj, lonobj : array-like
            1-D latitude / longitude arrays.
        dem_range : tuple of (float, float), optional
            (min, max) DEM elevation.  When provided, only the height
            layers spanning this range (plus padding) are loaded,
            reducing memory and computation.
        """
        self.offset = None
        self.interp = []
        self.latobj = latobj[:]
        self.lonobj = lonobj[:]

        hgts_full = np.asarray(hgtobj[:])
        data_full = np.asarray(inobj[:])

        if dem_range is not None:
            # Cubic vertical interp needs 2 extra layers on each side
            idx = _get_height_subset_indices(
                hgts_full, dem_range[0], dem_range[1], pad=2)
            self.hgts = hgts_full[idx]
            self.data = data_full[idx]
        else:
            self.hgts = hgts_full
            self.data = data_full

        self.createInterp()

    def createInterp(self):
        """Create interpolators."""
        self.offset = np.mean(self.data)
        for i in range(len(self.hgts)):
            self.interp.append(scipy.interpolate.RectBivariateSpline(
                self.latobj, self.lonobj, self.data[i] - self.offset))

    def __call__(self, line, pix, h):
        """Interpolate at a single point."""
        vals = np.array([x(line, pix)[0, 0] for x in self.interp])
        est = scipy.interpolate.interp1d(self.hgts, vals, kind='cubic')
        return est(h) + self.offset
