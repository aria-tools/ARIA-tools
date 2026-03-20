# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import scipy.interpolate


def _get_height_subset_indices(heights, dem_min, dem_max, pad=1):
    """Determine which height layers span the DEM elevation range.

    Parameters
    ----------
    heights : np.ndarray
        1-D array of height levels (ascending or descending).
    dem_min, dem_max : float
        Minimum / maximum DEM elevation.
    pad : int, optional
        Extra layers on each side of the bracket (default 1).

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
