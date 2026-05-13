"""
Mask caching utilities for ARIA-tools.

Provides a singleton cache for mask arrays to avoid repeatedly opening
and reading the same mask file during extraction workflows.

Author: ARIA-tools team
Copyright (c) 2026, California Institute of Technology.
"""

import logging
import threading
import numpy as np
from osgeo import gdal

LOGGER = logging.getLogger(__name__)


class MaskCache:
    """
    Thread-safe singleton cache for mask arrays.

    Masks are read once from disk and cached in memory for reuse across
    multiple products/layers in the same extraction run.

    Usage:
        # Get mask array (reads once, caches for reuse)
        mask_array = MaskCache.get('/path/to/mask.tif')

        # Clear cache between runs
        MaskCache.clear()

    Notes:
        - Thread-safe for parallel extraction
        - GDAL file handles are properly closed after reading
        - Returns None for None input (no-mask case)
        - Cache persists for lifetime of process unless cleared
    """

    _cache = {}
    _lock = threading.RLock()
    _stats = {'hits': 0, 'misses': 0}

    @classmethod
    def get(cls, maskfile):
        """
        Get mask array, reading from file only once.

        Parameters
        ----------
        maskfile : str or None
            Path to mask file. If None, returns None.

        Returns
        -------
        np.ndarray or None
            Mask array, or None if maskfile is None.

        Notes
        -----
        Thread-safe: Multiple threads can call simultaneously without
        duplicate reads or race conditions.
        """
        if maskfile is None:
            return None

        # Fast path: check cache without lock (common case)
        if maskfile in cls._cache:
            with cls._lock:
                cls._stats['hits'] += 1
            return cls._cache[maskfile]

        # Slow path: read from disk with lock
        with cls._lock:
            # Double-check after acquiring lock (another thread may have loaded)
            if maskfile in cls._cache:
                cls._stats['hits'] += 1
                return cls._cache[maskfile]

            # Read mask from disk
            LOGGER.debug('Loading mask from file: %s', maskfile)
            try:
                mask_ds = gdal.Open(maskfile, gdal.GA_ReadOnly)
                if mask_ds is None:
                    LOGGER.warning('Failed to open mask file: %s', maskfile)
                    return None

                # Read array
                mask_array = mask_ds.ReadAsArray()

                # Explicitly close GDAL dataset
                mask_ds = None

                # Cache the array
                cls._cache[maskfile] = mask_array
                cls._stats['misses'] += 1

                LOGGER.debug('Cached mask array: %s (shape: %s)',
                           maskfile, mask_array.shape)

                return mask_array

            except Exception as e:
                LOGGER.error('Error reading mask file %s: %s', maskfile, e)
                return None

    @classmethod
    def clear(cls):
        """
        Clear the mask cache.

        Call this between extraction runs to free memory and ensure
        fresh reads if mask files have changed.
        """
        with cls._lock:
            num_cached = len(cls._cache)
            cls._cache.clear()
            cls._stats['hits'] = 0
            cls._stats['misses'] = 0

            if num_cached > 0:
                LOGGER.debug('Cleared mask cache (%d masks freed)', num_cached)

    @classmethod
    def get_stats(cls):
        """
        Get cache statistics.

        Returns
        -------
        dict
            Dictionary with keys:
            - 'cached': Number of masks currently cached
            - 'hits': Total cache hits
            - 'misses': Total cache misses
            - 'hit_rate': Cache hit rate (0.0 to 1.0)
        """
        with cls._lock:
            total = cls._stats['hits'] + cls._stats['misses']
            hit_rate = cls._stats['hits'] / total if total > 0 else 0.0

            return {
                'cached': len(cls._cache),
                'hits': cls._stats['hits'],
                'misses': cls._stats['misses'],
                'hit_rate': hit_rate,
            }

    @classmethod
    def get_cached_files(cls):
        """
        Get list of currently cached mask files.

        Returns
        -------
        list of str
            Paths to cached mask files.
        """
        with cls._lock:
            return list(cls._cache.keys())


# Legacy function for backward compatibility
def get_mask_array(maskfile):
    """
    Get mask array using cache.

    Legacy wrapper around MaskCache.get() for backward compatibility.

    Parameters
    ----------
    maskfile : str or None
        Path to mask file.

    Returns
    -------
    np.ndarray or None
        Mask array.
    """
    return MaskCache.get(maskfile)
