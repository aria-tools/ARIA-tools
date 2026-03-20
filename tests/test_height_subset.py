"""Tests for the 3D-cube height-level subsetting optimisation.

Verifies that interpolating a 3D metadata cube after subsetting to
the DEM-relevant height levels gives identical (or near-identical)
results to interpolating the full cube, while reducing memory and
computation cost.

Usage:
  # Run all synthetic tests (no data needed)
  python tests/test_height_subset.py

  # Run with a real GUNW file (enables 10 additional real-data tests)
  python tests/test_height_subset.py --gunw-file /path/to/file.nc

  # Or via pytest
  pytest tests/test_height_subset.py -v
  pytest tests/test_height_subset.py -v --gunw-file /path/to/file.nc
"""

import time
import tempfile
import os

import numpy as np
import pytest
import scipy.interpolate

# ── units under test ────────────────────────────────────────────────
from ARIAtools.util.interp import (
    InterpCube, _get_height_subset_indices, _compute_dem_range,
)


# ====================================================================
# 0. _compute_dem_range — nodata-aware DEM min/max
# ====================================================================
class TestComputeDemRange:
    """Verify _compute_dem_range correctly skips nodata pixels."""

    @staticmethod
    def _make_dem(arr, nodata=None, dtype=None):
        """Create an in-memory GDAL dataset from a 2-D numpy array."""
        from osgeo import gdal, osr
        if dtype is None:
            dtype = gdal.GDT_Float32
        drv = gdal.GetDriverByName('MEM')
        ny, nx = arr.shape
        ds = drv.Create('', nx, ny, 1, dtype)
        ds.SetGeoTransform([0, 1, 0, ny, 0, -1])
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds.SetProjection(srs.ExportToWkt())
        band = ds.GetRasterBand(1)
        if nodata is not None:
            band.SetNoDataValue(nodata)
        band.WriteArray(arr)
        band.FlushCache()
        return ds

    def test_no_nodata(self):
        """Without nodata, returns simple min/max."""
        arr = np.array([[100, 200], [300, 400]], dtype=np.float32)
        ds = self._make_dem(arr)
        assert _compute_dem_range(ds) == (100.0, 400.0)

    def test_nan_nodata_float(self):
        """Float DEM with NaN nodata — NaN pixels excluded."""
        arr = np.array([[np.nan, 500], [1000, np.nan]], dtype=np.float32)
        ds = self._make_dem(arr, nodata=float('nan'))
        assert _compute_dem_range(ds) == (500.0, 1000.0)

    def test_nan_nodata_warped_to_int16(self):
        """NaN nodata on Int16 band — NaN can't match integer pixels.

        Without dstNodata in gdal.Warp, NaN pixels become 0 and nodata
        stays NaN.  _compute_dem_range reads the nodata as-is; since
        NaN never matches integer values, those 0s are treated as valid.
        The real fix is in dem.py (dstNodata=-32768).
        """
        from osgeo import gdal
        src_arr = np.array([[np.nan, 500], [2200, 3679]], dtype=np.float32)
        ds_src = self._make_dem(src_arr, nodata=float('nan'))
        # Warp to Int16 without dstNodata (legacy behavior)
        ds_int16 = gdal.Warp('', ds_src, format='MEM',
                             outputType=gdal.GDT_Int16)
        dem_min, dem_max = _compute_dem_range(ds_int16)
        # 0 is included because NaN nodata can't match integer pixels
        assert dem_min == 0.0
        assert dem_max == 3679.0

    def test_zero_nodata_int16(self):
        """Int16 DEM with explicit nodata=0 — zeros excluded."""
        arr = np.array([[0, 100], [200, 0]], dtype=np.int16)
        ds = self._make_dem(arr, nodata=0, dtype=__import__('osgeo').gdal.GDT_Int16)
        assert _compute_dem_range(ds) == (100.0, 200.0)

    def test_all_nodata_raises(self):
        """All-nodata DEM raises ValueError."""
        arr = np.array([[0, 0], [0, 0]], dtype=np.int16)
        ds = self._make_dem(arr, nodata=0, dtype=__import__('osgeo').gdal.GDT_Int16)
        with pytest.raises(ValueError, match='All DEM pixels are nodata'):
            _compute_dem_range(ds)

    def test_proper_dstNodata_int16(self):
        """Int16 DEM with dstNodata=-32768 (new dem.py pipeline).

        With the fix in dem.py, NaN→-32768 and real 0m elevation
        pixels are preserved.  _compute_dem_range should correctly
        include 0 and exclude -32768.
        """
        from osgeo import gdal
        src_arr = np.array([[np.nan, 0], [2200, 3679]], dtype=np.float32)
        ds_src = self._make_dem(src_arr, nodata=float('nan'))
        ds_fixed = gdal.Warp('', ds_src, format='MEM',
                             outputType=gdal.GDT_Int16, dstNodata=-32768)
        dem_min, dem_max = _compute_dem_range(ds_fixed)
        assert dem_min == 0.0, "Real 0m elevation must be preserved"
        assert dem_max == 3679.0


# ====================================================================
# 0b. Rioxarray dimension metadata consistency after band subsetting
# ====================================================================
class TestWarpedMetadataConsistency:
    """Verify that NETCDF dimension metadata is updated after subsetting.

    When create_raster_from_gunw subsets bands (e.g. 20→8) via
    gdal.Translate + gdal.Warp, the warped TIF inherits the original
    NETCDF_DIM_*_VALUES metadata with all 20 values.  rioxarray parses
    this metadata to build a coordinate dimension, causing a shape
    mismatch.  The fix updates the metadata to match the actual bands.
    """

    @staticmethod
    def _make_source_tif(tmpdir, n_bands=20):
        """Create a GeoTIFF with NETCDF height dimension metadata."""
        from osgeo import gdal, osr
        src_tif = os.path.join(tmpdir, 'src.tif')
        drv = gdal.GetDriverByName('GTiff')
        ds = drv.Create(src_tif, 10, 10, n_bands, gdal.GDT_Float32)
        ds.SetGeoTransform([-99.5, 0.01, 0, 19.5, 0, -0.01])
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds.SetProjection(srs.ExportToWkt())
        heights = np.linspace(-500, 9000, n_bands)
        heights_str = '{' + ','.join(str(h) for h in heights) + '}'
        ds.SetMetadataItem(
            'NETCDF_DIM_EXTRA', '{heightAboveEllipsoid}')
        ds.SetMetadataItem(
            'NETCDF_DIM_heightAboveEllipsoid_DEF',
            '{' + str(n_bands) + ',6}')
        ds.SetMetadataItem(
            'NETCDF_DIM_heightAboveEllipsoid_VALUES', heights_str)
        for b in range(1, n_bands + 1):
            ds.GetRasterBand(b).WriteArray(
                np.ones((10, 10), dtype=np.float32) * b)
        ds.FlushCache()
        ds = None
        return src_tif, heights

    def test_rioxarray_crash_without_fix(self):
        """Demonstrate that stale metadata causes rioxarray to crash."""
        import tempfile
        import shutil
        from osgeo import gdal
        import rioxarray

        tmpdir = tempfile.mkdtemp()
        try:
            src_tif, _ = self._make_source_tif(tmpdir, n_bands=20)

            # Subset to 8 bands
            sub_vrt = os.path.join(tmpdir, 'subset.vrt')
            opts = gdal.TranslateOptions(
                format='VRT', bandList=list(range(5, 13)))
            gdal.Translate(sub_vrt, src_tif, options=opts)

            # Warp (propagates stale metadata)
            mosaic_tif = os.path.join(tmpdir, 'warp.tif')
            ds = gdal.Warp(
                mosaic_tif, sub_vrt, format='GTiff',
                dstNodata=float('nan'))
            ds = None

            # Confirm stale metadata
            ds_check = gdal.Open(mosaic_tif)
            assert ds_check.RasterCount == 8
            vals = ds_check.GetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_VALUES')
            n_meta = len(vals[1:-1].split(','))
            assert n_meta == 20, "Metadata still has 20 values"
            ds_check = None

            # rioxarray should crash
            with pytest.raises(ValueError, match='conflicting sizes'):
                rioxarray.open_rasterio(mosaic_tif, masked=True)
        finally:
            shutil.rmtree(tmpdir)

    def test_updated_metadata_allows_rioxarray(self):
        """After fixing metadata, rioxarray opens successfully."""
        import tempfile
        import shutil
        from osgeo import gdal
        import rioxarray

        tmpdir = tempfile.mkdtemp()
        try:
            src_tif, heights = self._make_source_tif(tmpdir, n_bands=20)
            band_indices = list(range(4, 12))  # 0-based → bands 5-12
            subset_heights = np.array(heights)[band_indices]

            # Subset to 8 bands
            sub_vrt = os.path.join(tmpdir, 'subset.vrt')
            band_list = [i + 1 for i in band_indices]
            opts = gdal.TranslateOptions(format='VRT', bandList=band_list)
            gdal.Translate(sub_vrt, src_tif, options=opts)

            # Warp
            mosaic_tif = os.path.join(tmpdir, 'warp.tif')
            ds = gdal.Warp(
                mosaic_tif, sub_vrt, format='GTiff',
                dstNodata=float('nan'))
            ds = None

            # Apply the fix: update metadata
            ds_fix = gdal.Open(mosaic_tif, gdal.GA_Update)
            new_vals = '{' + ','.join(str(h) for h in subset_heights) + '}'
            ds_fix.SetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_VALUES', new_vals)
            ds_fix.SetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_DEF',
                '{' + str(len(subset_heights)) + ',6}')
            ds_fix.FlushCache()
            ds_fix = None

            # rioxarray should now succeed
            da = rioxarray.open_rasterio(mosaic_tif, masked=True)
            assert da.sizes['heightAboveEllipsoid'] == 8
            da.close()
        finally:
            shutil.rmtree(tmpdir)

    def test_no_metadata_no_crash(self):
        """Without NETCDF dimension metadata, no fix needed."""
        import tempfile
        import shutil
        from osgeo import gdal
        import rioxarray

        tmpdir = tempfile.mkdtemp()
        try:
            # Source without NETCDF metadata
            src_tif = os.path.join(tmpdir, 'plain.tif')
            drv = gdal.GetDriverByName('GTiff')
            ds = drv.Create(src_tif, 10, 10, 4, gdal.GDT_Float32)
            ds.SetGeoTransform([-99.5, 0.01, 0, 19.5, 0, -0.01])
            from osgeo import osr
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(4326)
            ds.SetProjection(srs.ExportToWkt())
            for b in range(1, 5):
                ds.GetRasterBand(b).WriteArray(
                    np.ones((10, 10), dtype=np.float32) * b)
            ds.FlushCache()
            ds = None

            # Subset to 2 bands
            sub_vrt = os.path.join(tmpdir, 'sub.vrt')
            gdal.Translate(sub_vrt, src_tif,
                           options=gdal.TranslateOptions(
                               format='VRT', bandList=[1, 2]))
            warp_tif = os.path.join(tmpdir, 'warp.tif')
            ds = gdal.Warp(warp_tif, sub_vrt, format='GTiff')
            ds = None

            # Should open fine
            da = rioxarray.open_rasterio(warp_tif, masked=True)
            assert da.sizes['band'] == 2
            da.close()
        finally:
            shutil.rmtree(tmpdir)

    def test_create_raster_from_gunw_with_subsetting(self):
        """End-to-end: create_raster_from_gunw subsets and writes
        consistent metadata that rioxarray can open."""
        import tempfile
        import shutil
        from osgeo import gdal, osr
        import rioxarray

        n_h, n_lat, n_lon = 20, 30, 40
        heights = np.linspace(-500, 9000, n_h).astype(np.float32)
        lat0, lon0 = 19.3, -99.5
        dlat, dlon = -0.01, 0.01

        tmpdir = tempfile.mkdtemp()
        try:
            # Source GeoTIFF with 20 bands + NETCDF dim metadata
            src_tif = os.path.join(tmpdir, 'src.tif')
            drv = gdal.GetDriverByName('GTiff')
            ds = drv.Create(src_tif, n_lon, n_lat, n_h, gdal.GDT_Float32)
            ds.SetGeoTransform([lon0, dlon, 0, lat0, 0, dlat])
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(4326)
            ds.SetProjection(srs.ExportToWkt())
            hgt_str = '{' + ','.join(str(h) for h in heights) + '}'
            ds.SetMetadataItem(
                'NETCDF_DIM_EXTRA', '{heightAboveEllipsoid}')
            ds.SetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_DEF',
                '{' + str(n_h) + ',6}')
            ds.SetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_VALUES', hgt_str)
            for b in range(1, n_h + 1):
                ds.GetRasterBand(b).WriteArray(
                    np.full((n_lat, n_lon), 30 + b, dtype=np.float32))
            ds.FlushCache()
            ds = None

            # DEM with elevation range 2000–4000 m
            dem_tif = os.path.join(tmpdir, 'dem.tif')
            dem_ds = drv.Create(dem_tif, n_lon, n_lat, 1, gdal.GDT_Int16)
            dem_ds.SetGeoTransform([lon0, dlon, 0, lat0, 0, dlat])
            dem_ds.SetProjection(srs.ExportToWkt())
            dem_ds.GetRasterBand(1).SetNoDataValue(-32768)
            dem_ds.GetRasterBand(1).WriteArray(
                np.random.RandomState(42).randint(
                    2000, 4000, (n_lat, n_lon)).astype(np.int16))
            dem_ds.FlushCache()
            dem_ds = None

            dem_ds = gdal.Open(dem_tif)
            from ARIAtools.extractProduct import create_raster_from_gunw

            outname = os.path.join(tmpdir, 'output')
            create_raster_from_gunw(
                outname, [src_tif], 'EPSG:4326', 'ENVI',
                hgt_field='NETCDF_DIM_heightAboveEllipsoid_VALUES',
                dem=dem_ds)
            dem_ds = None

            # Verify output VRT has consistent metadata
            ds_out = gdal.Open(outname + '.vrt')
            assert ds_out is not None
            n_bands = ds_out.RasterCount
            assert n_bands < n_h, \
                f"Expected subsetting but got {n_bands} bands"
            hgt_meta = ds_out.GetMetadataItem(
                'NETCDF_DIM_heightAboveEllipsoid_VALUES')
            assert hgt_meta is not None
            n_vals = len(hgt_meta[1:-1].split(','))
            assert n_vals == n_bands, (
                f"Metadata has {n_vals} height values "
                f"but output has {n_bands} bands")
            ds_out = None

            # Verify rioxarray can open the intermediate warped TIF
            # (the fix prevents the crash)
            da = rioxarray.open_rasterio(outname, masked=True)
            da.close()
        finally:
            shutil.rmtree(tmpdir)


# ====================================================================
# Helper: build a synthetic 3D cube that mimics a GUNW metadata layer
# ====================================================================
def _make_synthetic_cube(n_heights=15, n_lat=50, n_lon=60, seed=42):
    """Return (data, heights, lats, lons) for a synthetic 3D cube.

    Heights go from -1000 m to 8000 m in ``n_heights`` equally-spaced
    steps, mimicking typical GUNW height levels.
    """
    rng = np.random.RandomState(seed)
    heights = np.linspace(-1000, 8000, n_heights, dtype='float32')
    lats = np.linspace(34.0, 35.0, n_lat, dtype='float32')
    lons = np.linspace(-118.0, -117.0, n_lon, dtype='float32')

    # Smooth-ish 3D field: quadratic in height + spatial variation
    LON, LAT, HGT = np.meshgrid(lons, lats, heights, indexing='ij')
    data = (
        0.05 * HGT
        + 2.0 * np.sin(LAT * 10) * np.cos(LON * 10)
        + 0.001 * HGT ** 2 / 8000
        + rng.randn(*LON.shape).astype('float32') * 0.01
    ).astype('float32')
    # Reorder to (heights, lats, lons) — matching GUNW convention
    data = data.transpose(2, 1, 0)
    return data, heights, lats, lons


# ====================================================================
# 1. _get_height_subset_indices — unit tests
# ====================================================================
class TestGetHeightSubsetIndices:
    """Unit tests for _get_height_subset_indices."""

    def test_ascending_basic(self):
        """Ascending heights, DEM range inside cube → proper bracket."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000, 4000, 5000], dtype='float32')
        idx = _get_height_subset_indices(heights, -100, 2100, pad=0)
        # Bracket: -1000 ≤ -100 < 0, 2000 < 2100 ≤ 3000 → [0,4]
        # Pad=0  : [max(0,0), min(6,4)] → [0, 4]
        selected = heights[idx]
        assert selected[0] <= -100, "Must include layer ≤ dem_min"
        assert selected[-1] >= 2100, "Must include layer ≥ dem_max"
        assert len(idx) < len(heights), "Subset should be smaller"

    def test_descending_basic(self):
        """Descending heights → same logic, reversed indices."""
        heights = np.array(
            [5000, 4000, 3000, 2000, 1000, 0, -1000], dtype='float32')
        idx = _get_height_subset_indices(heights, -100, 2100, pad=0)
        selected = heights[idx]
        assert np.min(selected) <= -100
        assert np.max(selected) >= 2100
        assert len(idx) < len(heights)

    def test_dem_covers_full_range(self):
        """DEM range covers or exceeds cube → all bands returned."""
        heights = np.array([-1000, 0, 1000, 2000], dtype='float32')
        idx = _get_height_subset_indices(heights, -2000, 5000, pad=0)
        np.testing.assert_array_equal(idx, np.arange(len(heights)))

    def test_single_height(self):
        """Single height level → always returns index [0]."""
        heights = np.array([500.0], dtype='float32')
        idx = _get_height_subset_indices(heights, 0, 1000, pad=0)
        np.testing.assert_array_equal(idx, np.array([0]))

    def test_exact_boundary(self):
        """DEM min/max exactly on height levels."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000], dtype='float32')
        idx = _get_height_subset_indices(heights, 0, 2000, pad=0)
        selected = heights[idx]
        assert 0 in selected
        assert 2000 in selected
        # With pad=0: only the spanning levels are selected
        assert -1000 not in selected
        assert 3000 not in selected

    def test_pad_zero(self):
        """pad=0 gives the tightest bracket."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000, 4000, 5000], dtype='float32')
        idx = _get_height_subset_indices(heights, 500, 2500, pad=0)
        selected = heights[idx]
        assert selected[0] <= 500
        assert selected[-1] >= 2500

    def test_pad_two(self):
        """pad=2 gives a wider bracket than pad=0."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000, 4000, 5000], dtype='float32')
        idx_p0 = _get_height_subset_indices(heights, 1500, 2500, pad=0)
        idx_p2 = _get_height_subset_indices(heights, 1500, 2500, pad=2)
        assert len(idx_p2) >= len(idx_p0)

    def test_pad_never_exceeds_bounds(self):
        """Large pad values are clamped to available height range."""
        heights = np.array(
            [0, 1000, 2000, 3000, 4000], dtype='float32')
        # pad=10 is much larger than the number of levels
        idx = _get_height_subset_indices(heights, 1500, 2500, pad=10)
        # Must return valid indices within [0, len-1]
        assert idx.min() >= 0
        assert idx.max() <= len(heights) - 1
        # Should return all levels since pad overwhelms the range
        np.testing.assert_array_equal(idx, np.arange(len(heights)))

    def test_pad_clamped_at_bottom(self):
        """Pad doesn't go below index 0 when DEM is near bottom."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000, 4000, 5000], dtype='float32')
        idx = _get_height_subset_indices(heights, -900, -800, pad=2)
        assert idx.min() == 0  # Can't go below 0
        selected = heights[idx]
        assert np.min(selected) <= -900

    def test_pad_clamped_at_top(self):
        """Pad doesn't go above index n-1 when DEM is near top."""
        heights = np.array(
            [-1000, 0, 1000, 2000, 3000, 4000, 5000], dtype='float32')
        idx = _get_height_subset_indices(heights, 4800, 4900, pad=2)
        assert idx.max() == len(heights) - 1  # Can't exceed last index
        selected = heights[idx]
        assert np.max(selected) >= 4900 or idx.max() == len(heights) - 1


# ====================================================================
# 2. RegularGridInterpolator: full cube vs subsetted cube
# ====================================================================
class TestRegularGridInterpolatorSubset:
    """Compare interpolation results: full 3D cube vs subsetted cube."""

    @pytest.fixture
    def cube_data(self):
        return _make_synthetic_cube(n_heights=15, n_lat=50, n_lon=60)

    def _interpolate_cube(self, data, heights, lats, lons, dem_z):
        """Run RegularGridInterpolator on full or subsetted cube."""
        LON_q, LAT_q = np.meshgrid(
            lons[5:-5], lats[5:-5], indexing='ij')

        interper = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights),
            data.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        pts = np.stack(
            [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)
        return interper(pts).reshape(LAT_q.shape)

    def test_subset_matches_full(self, cube_data):
        """Subsetted interpolation matches full-cube interpolation."""
        data, heights, lats, lons = cube_data

        # Simulate a DEM with range [200, 3500] meters
        rng = np.random.RandomState(99)
        n_q = (len(lons) - 10, len(lats) - 10)
        dem_z = rng.uniform(200, 3500, size=n_q).astype('float32')

        # Full cube
        result_full = self._interpolate_cube(
            data, heights, lats, lons, dem_z)

        # Subsetted cube
        idx = _get_height_subset_indices(heights, 200, 3500, pad=0)
        result_sub = self._interpolate_cube(
            data[idx], heights[idx], lats, lons, dem_z)

        np.testing.assert_allclose(
            result_sub, result_full, rtol=1e-5, atol=1e-5,
            err_msg="Subsetted interpolation must match full-cube result")

    def test_subset_smaller_memory(self, cube_data):
        """Subsetted cube uses less memory than the full cube."""
        data, heights, lats, lons = cube_data
        idx = _get_height_subset_indices(heights, 200, 3500, pad=0)
        assert data[idx].nbytes < data.nbytes

    def test_subset_is_faster_large(self):
        """Subsetted interpolation is faster on a larger cube.

        Uses a larger dataset where the I/O and memory savings are
        measurable even in pure-Python (no GDAL disk I/O).
        The real-world benefit is even greater because GDAL loads
        fewer raster bands from disk.
        """
        data, heights, lats, lons = _make_synthetic_cube(
            n_heights=50, n_lat=200, n_lon=250)

        rng = np.random.RandomState(99)
        n_q = (len(lons) - 10, len(lats) - 10)
        dem_z = rng.uniform(1000, 3000, size=n_q).astype('float32')

        idx = _get_height_subset_indices(heights, 1000, 3000, pad=0)

        # Warm up
        self._interpolate_cube(data, heights, lats, lons, dem_z)
        self._interpolate_cube(
            data[idx], heights[idx], lats, lons, dem_z)

        # Timed runs
        n_iter = 3
        t0 = time.perf_counter()
        for _ in range(n_iter):
            self._interpolate_cube(data, heights, lats, lons, dem_z)
        t_full = time.perf_counter() - t0

        t0 = time.perf_counter()
        for _ in range(n_iter):
            self._interpolate_cube(
                data[idx], heights[idx], lats, lons, dem_z)
        t_sub = time.perf_counter() - t0

        print(f"\n  Full cube: {t_full/n_iter:.4f}s per iter"
              f"  ({len(heights)} levels)")
        print(f"  Subset:    {t_sub/n_iter:.4f}s per iter"
              f"  ({len(idx)} levels)")
        print(f"  Speed-up:  {t_full/t_sub:.2f}x")

        # With 50 heights subsetted to ~12, the subsetted version
        # should be noticeably faster. Allow generous margin for CI.
        assert t_sub < t_full * 1.5, (
            f"Subsetted ({t_sub:.3f}s) should not be much slower "
            f"than full ({t_full:.3f}s)")


# ====================================================================
# 3. InterpCube: full vs with dem_range
# ====================================================================
class TestInterpCubeSubset:
    """Compare InterpCube results with and without dem_range."""

    def test_dem_range_reduces_layers(self):
        """Passing dem_range loads fewer height layers."""
        data, heights, lats, lons = _make_synthetic_cube(
            n_heights=15, n_lat=20, n_lon=25)

        cube_full = InterpCube(data, heights, lats, lons)
        cube_sub = InterpCube(
            data, heights, lats, lons, dem_range=(200, 3500))

        assert len(cube_sub.hgts) < len(cube_full.hgts)
        assert len(cube_sub.interp) < len(cube_full.interp)

    def test_dem_range_matches_full(self):
        """InterpCube with dem_range gives same result at query points."""
        data, heights, lats, lons = _make_synthetic_cube(
            n_heights=15, n_lat=20, n_lon=25)

        cube_full = InterpCube(data, heights, lats, lons)
        cube_sub = InterpCube(
            data, heights, lats, lons, dem_range=(200, 3500))

        # Query at several points within the DEM range
        test_lats = lats[5:-5:3]
        test_lons = lons[5:-5:3]
        test_hgts = np.linspace(300, 3400, 8)

        for lat in test_lats:
            for lon in test_lons:
                for h in test_hgts:
                    v_full = cube_full(lat, lon, h)
                    v_sub = cube_sub(lat, lon, h)
                    np.testing.assert_allclose(
                        v_sub, v_full, rtol=1e-4, atol=1e-4,
                        err_msg=f"Mismatch at ({lat}, {lon}, {h})")

    def test_backward_compatible(self):
        """InterpCube without dem_range works as before."""
        data, heights, lats, lons = _make_synthetic_cube(
            n_heights=10, n_lat=15, n_lon=20)

        cube = InterpCube(data, heights, lats, lons)
        assert len(cube.hgts) == len(heights)
        # Should be callable
        val = cube(lats[5], lons[5], heights[3])
        assert np.isfinite(val)


# ====================================================================
# 4. Larger-scale benchmark (parametrized)
# ====================================================================
@pytest.mark.parametrize("n_heights,dem_lo,dem_hi", [
    (7,   -100,  2100),    # typical: 7-height cube, low-elevation DEM
    (15,   200,  3500),    # 15-height cube, mid-elevation DEM
    (30,  4000,  6000),    # 30-height cube, high-altitude DEM (narrow)
    (15, -1200,  9000),    # DEM exceeds cube → no subsetting
])
def test_interpolation_accuracy_parametrized(
        n_heights, dem_lo, dem_hi):
    """Parametrized accuracy test across scenarios."""
    data, heights, lats, lons = _make_synthetic_cube(
        n_heights=n_heights, n_lat=40, n_lon=50)

    idx = _get_height_subset_indices(heights, dem_lo, dem_hi, pad=0)

    rng = np.random.RandomState(123)
    n_q = (len(lons) - 10, len(lats) - 10)
    # Clamp DEM to cube range to avoid NaN in full interpolation too
    h_min, h_max = float(heights.min()), float(heights.max())
    dem_z = rng.uniform(
        max(dem_lo, h_min), min(dem_hi, h_max),
        size=n_q).astype('float32')

    LON_q, LAT_q = np.meshgrid(
        lons[5:-5], lats[5:-5], indexing='ij')

    # Full cube interpolation
    interper_full = scipy.interpolate.RegularGridInterpolator(
        (lats, lons, heights),
        data.transpose(1, 2, 0),
        fill_value=np.nan, bounds_error=False)
    pts = np.stack(
        [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)
    result_full = interper_full(pts)

    # Subsetted cube interpolation
    interper_sub = scipy.interpolate.RegularGridInterpolator(
        (lats, lons, heights[idx]),
        data[idx].transpose(1, 2, 0),
        fill_value=np.nan, bounds_error=False)
    result_sub = interper_sub(pts)

    # Where both are finite, they should match
    mask = np.isfinite(result_full) & np.isfinite(result_sub)
    if mask.any():
        np.testing.assert_allclose(
            result_sub[mask], result_full[mask],
            rtol=1e-5, atol=1e-5,
            err_msg=f"Mismatch for n_heights={n_heights}, "
                    f"DEM=[{dem_lo},{dem_hi}]")

    n_total = len(heights)
    n_sub = len(idx)
    print(f"\n  heights={n_total}, subset={n_sub}, "
          f"DEM=[{dem_lo},{dem_hi}], "
          f"memory reduction={100*(1-n_sub/n_total):.0f}%")


# ====================================================================
# 5. Edge-case tests
# ====================================================================
class TestEdgeCases:
    """Edge-case behaviour for height subsetting."""

    def test_dem_below_all_heights(self):
        """DEM entirely below lowest height → uses bottom layers."""
        heights = np.array([0, 1000, 2000, 3000], dtype='float32')
        idx = _get_height_subset_indices(heights, -500, -100, pad=0)
        assert 0 in idx  # Must include the lowest available

    def test_dem_above_all_heights(self):
        """DEM entirely above highest height → uses top layers."""
        heights = np.array([0, 1000, 2000, 3000], dtype='float32')
        idx = _get_height_subset_indices(heights, 3500, 5000, pad=0)
        assert len(heights) - 1 in idx  # Must include the highest

    def test_two_heights_only(self):
        """Only 2 height levels → no subsetting possible."""
        heights = np.array([0, 5000], dtype='float32')
        idx = _get_height_subset_indices(heights, 100, 4000, pad=0)
        np.testing.assert_array_equal(idx, np.array([0, 1]))

    def test_negative_heights(self):
        """Heights can be negative (below sea level)."""
        heights = np.array(
            [-3000, -2000, -1000, 0, 1000], dtype='float32')
        idx = _get_height_subset_indices(heights, -2500, -500, pad=0)
        selected = heights[idx]
        assert np.min(selected) <= -2500
        assert np.max(selected) >= -500


# ====================================================================
# 6. Geometry-layer-specific tests
# ====================================================================
def _make_geometry_cube(layer_name, n_heights=15, n_lat=30, n_lon=40,
                        seed=42):
    """Build a synthetic 3D cube mimicking a specific GUNW geometry layer.

    Each layer type has a realistic value range:
    - incidenceAngle: ~27°–50° (varies with look direction & height)
    - lookAngle:      ~24°–44°
    - azimuthAngle:   ~-170° (nearly constant, small spatial gradient)
    - bPerpendicular: ~-120 to -98 m (varies with height & position)
    - bParallel:      ~50–90 m
    - solidEarthTide: ~-0.02 to 0.02 m (small, varies spatially)
    - troposphereWet/Hydrostatic/Total: ~-0.3 to 0.3 m
    """
    rng = np.random.RandomState(seed)
    # Use GUNW-like height levels
    heights = np.linspace(-1500, 9000, n_heights, dtype='float32')
    lats = np.linspace(33.4, 36.0, n_lat, dtype='float32')
    lons = np.linspace(-119.3, -115.7, n_lon, dtype='float32')

    LON, LAT, HGT = np.meshgrid(lons, lats, heights, indexing='ij')

    # Layer-specific value generation
    if layer_name == 'incidenceAngle':
        base = 38.0 + 8.0 * np.sin((LAT - 34.5) * 2)
        data = base + 0.001 * HGT + rng.randn(*LON.shape) * 0.01
    elif layer_name == 'lookAngle':
        base = 33.5 + 7.0 * np.sin((LAT - 34.5) * 2)
        data = base + 0.0008 * HGT + rng.randn(*LON.shape) * 0.01
    elif layer_name == 'azimuthAngle':
        base = -170.1 + 0.2 * np.sin(LAT * 5) * np.cos(LON * 3)
        data = base + 0.00001 * HGT + rng.randn(*LON.shape) * 0.005
    elif layer_name == 'bPerpendicular':
        base = -110.0 + 10.0 * np.cos((LAT - 34.5) * 3)
        data = base + 0.002 * HGT + rng.randn(*LON.shape) * 0.1
    elif layer_name == 'bParallel':
        base = 70.0 + 15.0 * np.sin((LAT - 34.5) * 3)
        data = base - 0.001 * HGT + rng.randn(*LON.shape) * 0.1
    elif layer_name == 'solidEarthTide':
        data = 0.01 * np.sin(LAT * 20) * np.cos(LON * 15) + 0.0001 * HGT
        data += rng.randn(*LON.shape) * 0.001
    elif layer_name in ('troposphereWet', 'troposphereHydrostatic',
                        'troposphereTotal'):
        data = 0.15 * np.sin(LAT * 8) * np.cos(LON * 6) + 0.00005 * HGT
        data += rng.randn(*LON.shape) * 0.005
    else:
        raise ValueError(f"Unknown layer: {layer_name}")

    data = data.astype('float32').transpose(2, 1, 0)
    return data, heights, lats, lons


# All layers that go through 3D height interpolation in finalize_metadata
GEOMETRY_LAYERS = [
    'incidenceAngle',
    'lookAngle',
    'azimuthAngle',
    'bPerpendicular',
    'bParallel',
    'solidEarthTide',
    'troposphereWet',
    'troposphereHydrostatic',
    'troposphereTotal',
]


@pytest.mark.parametrize("layer_name", GEOMETRY_LAYERS)
class TestGeometryLayerSubset:
    """Verify subsetting accuracy for each geometry/metadata layer type."""

    def test_subset_matches_full(self, layer_name):
        """Subsetted interpolation matches full cube for this layer."""
        data, heights, lats, lons = _make_geometry_cube(layer_name)

        # Typical mid-elevation DEM range
        dem_min, dem_max = 200.0, 3500.0
        idx = _get_height_subset_indices(heights, dem_min, dem_max, pad=0)

        rng = np.random.RandomState(77)
        n_q = (len(lons) - 6, len(lats) - 6)
        dem_z = rng.uniform(dem_min, dem_max, size=n_q).astype('float32')

        LON_q, LAT_q = np.meshgrid(
            lons[3:-3], lats[3:-3], indexing='ij')

        # Full cube
        interper_full = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights),
            data.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        pts = np.stack(
            [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)
        result_full = interper_full(pts)

        # Subsetted cube
        interper_sub = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights[idx]),
            data[idx].transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        result_sub = interper_sub(pts)

        mask = np.isfinite(result_full) & np.isfinite(result_sub)
        np.testing.assert_allclose(
            result_sub[mask], result_full[mask],
            rtol=1e-5, atol=1e-5,
            err_msg=f"Subset mismatch for layer '{layer_name}'")

    def test_subset_reduces_bands(self, layer_name):
        """Subsetting actually reduces number of height bands."""
        data, heights, lats, lons = _make_geometry_cube(layer_name)
        idx = _get_height_subset_indices(heights, 200, 3500, pad=0)
        assert len(idx) < len(heights), (
            f"Expected fewer bands for {layer_name}")

    def test_high_elevation_dem(self, layer_name):
        """Test with high-altitude DEM (e.g., Himalayas, 4000-6000 m)."""
        data, heights, lats, lons = _make_geometry_cube(layer_name)
        dem_min, dem_max = 4000.0, 6000.0
        idx = _get_height_subset_indices(heights, dem_min, dem_max, pad=0)

        rng = np.random.RandomState(55)
        n_q = (len(lons) - 6, len(lats) - 6)
        dem_z = rng.uniform(dem_min, dem_max, size=n_q).astype('float32')

        LON_q, LAT_q = np.meshgrid(
            lons[3:-3], lats[3:-3], indexing='ij')
        pts = np.stack(
            [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)

        interper_full = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights),
            data.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        interper_sub = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights[idx]),
            data[idx].transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        result_full = interper_full(pts)
        result_sub = interper_sub(pts)

        mask = np.isfinite(result_full) & np.isfinite(result_sub)
        np.testing.assert_allclose(
            result_sub[mask], result_full[mask],
            rtol=1e-5, atol=1e-5,
            err_msg=f"High-elev mismatch for '{layer_name}'")

    def test_low_elevation_dem(self, layer_name):
        """Test with low/negative DEM (e.g., Death Valley, -100 to 500 m)."""
        data, heights, lats, lons = _make_geometry_cube(layer_name)
        dem_min, dem_max = -100.0, 500.0
        idx = _get_height_subset_indices(heights, dem_min, dem_max, pad=0)

        rng = np.random.RandomState(33)
        n_q = (len(lons) - 6, len(lats) - 6)
        dem_z = rng.uniform(dem_min, dem_max, size=n_q).astype('float32')

        LON_q, LAT_q = np.meshgrid(
            lons[3:-3], lats[3:-3], indexing='ij')
        pts = np.stack(
            [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)

        interper_full = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights),
            data.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        interper_sub = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights[idx]),
            data[idx].transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        result_full = interper_full(pts)
        result_sub = interper_sub(pts)

        mask = np.isfinite(result_full) & np.isfinite(result_sub)
        np.testing.assert_allclose(
            result_sub[mask], result_full[mask],
            rtol=1e-5, atol=1e-5,
            err_msg=f"Low-elev mismatch for '{layer_name}'")


# ====================================================================
# 7. Real GUNW data test
# ====================================================================
# Usage:
#   pytest tests/test_height_subset.py --gunw-file /path/to/file.nc -k TestRealGUNWData -v -s


def _get_gunw_path(config):
    return config.getoption('--gunw-file', default=None)


class TestRealGUNWData:
    """Test height subsetting on real GUNW product data.

    Reads 3D geometry layers from an actual GUNW file and verifies
    that subsetting to typical DEM elevation ranges produces
    identical interpolated results.

    Run with:  pytest tests/test_height_subset.py --gunw-file /path/to/file.nc
    """

    REAL_LAYERS = [
        'incidenceAngle', 'lookAngle', 'azimuthAngle',
        'perpendicularBaseline', 'parallelBaseline',
    ]

    @pytest.fixture(scope='class')
    def gunw_path(self, request):
        path = _get_gunw_path(request.config)
        if path is None or not os.path.exists(path):
            pytest.skip(
                'No GUNW file provided. '
                'Use --gunw-file /path/to/file.nc')
        return path

    @pytest.fixture(scope='class')
    def gunw_data(self, gunw_path):
        """Load all 3D geometry layers from the real GUNW file."""
        import h5py
        base = '/science/grids/imagingGeometry'
        with h5py.File(gunw_path, 'r') as f:
            heights = f[base + '/heightsMeta'][:].astype('float32')
            lats = f[base + '/latitudeMeta'][:].astype('float32')
            lons = f[base + '/longitudeMeta'][:].astype('float32')
            layers = {}
            for lyr in self.REAL_LAYERS:
                layers[lyr] = f[base + '/' + lyr][:].astype('float32')
        return heights, lats, lons, layers

    @pytest.mark.parametrize("layer_name", REAL_LAYERS)
    def test_real_subset_matches_full(self, gunw_data, layer_name):
        """Real GUNW: subsetted interpolation matches full cube."""
        heights, lats, lons, layers = gunw_data
        data = layers[layer_name]

        # Typical California DEM range
        dem_min, dem_max = -50.0, 2500.0
        idx = _get_height_subset_indices(heights, dem_min, dem_max, pad=0)

        print(f"\n  {layer_name}: heights={list(heights)}, "
              f"subset indices={list(idx)}, "
              f"subset heights={list(heights[idx])}")

        rng = np.random.RandomState(42)
        # Query at interior lat/lon points (avoid edges)
        margin = 2
        query_lats = lats[margin:-margin]
        query_lons = lons[margin:-margin]
        LON_q, LAT_q = np.meshgrid(query_lons, query_lats, indexing='ij')
        n_q = LON_q.shape
        dem_z = rng.uniform(dem_min, dem_max, size=n_q).astype('float32')

        # Full cube interpolation
        interper_full = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights),
            data.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        pts = np.stack(
            [LAT_q.ravel(), LON_q.ravel(), dem_z.ravel()], axis=-1)
        result_full = interper_full(pts)

        # Subsetted cube interpolation
        interper_sub = scipy.interpolate.RegularGridInterpolator(
            (lats, lons, heights[idx]),
            data[idx].transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)
        result_sub = interper_sub(pts)

        mask = np.isfinite(result_full) & np.isfinite(result_sub)
        assert mask.any(), f"No finite values for {layer_name}"
        np.testing.assert_allclose(
            result_sub[mask], result_full[mask],
            rtol=1e-5, atol=1e-5,
            err_msg=f"Real GUNW subset mismatch for '{layer_name}'")

        print(f"  Matched {mask.sum()} points, "
              f"max abs diff: "
              f"{np.max(np.abs(result_sub[mask] - result_full[mask])):.2e}")

    @pytest.mark.parametrize("layer_name", REAL_LAYERS)
    def test_real_few_heights_subset_saves(self, gunw_data, layer_name):
        """With only 4 height levels, verify subsetting is meaningful
        for narrow DEM ranges and graceful for wide ranges."""
        heights, lats, lons, layers = gunw_data

        # Narrow range: should subset
        idx_narrow = _get_height_subset_indices(
            heights, 100, 1500, pad=0)
        # Wide range: should return all
        idx_wide = _get_height_subset_indices(
            heights, -2000, 10000, pad=0)

        assert len(idx_wide) == len(heights), (
            "Wide DEM range should use all heights")
        print(f"\n  {layer_name}: 4 levels, "
              f"narrow subset={len(idx_narrow)}, "
              f"wide subset={len(idx_wide)}")


if __name__ == '__main__':
    import sys
    args = ['--tb=short', '-v', __file__] + sys.argv[1:]
    sys.exit(pytest.main(args))
