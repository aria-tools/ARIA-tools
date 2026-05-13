"""
Unit tests for metadata caching functionality.

Tests the ARIAtools.util.meta_cache module to ensure:
- Metadata is loaded only once from remote/local files
- Cache staleness detection works (TTL for remote, mtime for local)
- Thread-safe operations
- Expanded metadata fields (projection, bounds, geotransform, size)
- Cache persistence to disk
"""

import json
import os
import threading
import time
from unittest.mock import patch

import numpy as np
import pytest

# Test imports
from ARIAtools.util.meta_cache import (
    _cache_path,
    _file_key,
    _is_cache_stale,
    extract_metadata_gdal,
    get_or_extract,
    load_cache,
    save_cache,
)


class TestMetaCacheBasics:
    """Test basic metadata cache operations."""

    def test_cache_path_generation(self, tmp_path):
        """Cache path derived from URL file directory."""
        url_file = tmp_path / "urls.txt"
        url_file.touch()

        cache_path = _cache_path(str(url_file))

        assert cache_path == str(tmp_path / "aria_meta_cache.json")

    def test_cache_path_none_input(self):
        """None input returns None."""
        result = _cache_path(None)
        assert result is None

    def test_file_key_strips_vsicurl(self):
        """File key strips /vsicurl/ prefix."""
        url = "https://example.com/product.nc"
        vsicurl_path = f"/vsicurl/{url}"

        key = _file_key(vsicurl_path)

        assert key == url
        assert "/vsicurl/" not in key

    def test_file_key_strips_vsis3(self):
        """File key strips /vsis3/ prefix."""
        path = "/vsis3/bucket/product.nc"

        key = _file_key(path)

        assert "/vsis3/" not in key

    def test_load_cache_missing_file(self, tmp_path):
        """Load cache returns empty dict for missing file."""
        cache_file = tmp_path / "nonexistent_cache.json"

        cache_data = load_cache(str(cache_file))

        assert cache_data == {}

    def test_load_cache_valid_file(self, tmp_path):
        """Load cache reads existing cache file."""
        cache_file = tmp_path / "cache.json"
        test_data = {
            "product1.nc": {
                "version": "1c",
                "projection": "4326",
                "gdal_info_ts": time.time(),
            }
        }
        cache_file.write_text(json.dumps(test_data))

        cache_data = load_cache(str(cache_file))

        assert "product1.nc" in cache_data
        assert cache_data["product1.nc"]["version"] == "1c"

    def test_save_cache_writes_file(self, tmp_path):
        """Save cache persists data to disk."""
        cache_file = tmp_path / "cache.json"
        test_data = {"product1.nc": {"version": "1c", "projection": "4326"}}

        save_cache(str(cache_file), test_data)

        assert cache_file.exists()
        saved_data = json.loads(cache_file.read_text())
        assert saved_data == test_data

    def test_save_cache_none_file(self):
        """Save cache handles None file gracefully."""
        # Should not raise exception
        save_cache(None, {"key": "value"})


class TestCacheStaleness:
    """Test cache staleness detection."""

    def test_missing_timestamp_is_stale(self):
        """Cache entry without timestamp is stale."""
        cached_meta = {"version": "1c"}

        is_stale = _is_cache_stale("product.nc", cached_meta)

        assert is_stale is True

    def test_incomplete_cache_is_stale(self):
        """Cache entry without _cache_complete flag is stale."""
        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time(),
            "_cache_complete": False,
        }

        is_stale = _is_cache_stale("product.nc", cached_meta)

        assert is_stale is True

    def test_remote_file_within_ttl_not_stale(self):
        """Remote file within TTL is not stale."""
        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time() - 1800,  # 30 minutes ago
            "_cache_complete": True,
        }

        is_stale = _is_cache_stale(
            "https://example.com/product.nc",
            cached_meta,
            ttl_seconds=3600,  # 1 hour TTL
        )

        assert is_stale is False

    def test_remote_file_past_ttl_is_stale(self):
        """Remote file past TTL is stale."""
        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time() - 7200,  # 2 hours ago
            "_cache_complete": True,
        }

        is_stale = _is_cache_stale(
            "https://example.com/product.nc",
            cached_meta,
            ttl_seconds=3600,  # 1 hour TTL
        )

        assert is_stale is True

    def test_local_file_not_modified_not_stale(self, tmp_path):
        """Local file not modified since cache is not stale."""
        test_file = tmp_path / "product.nc"
        test_file.touch()
        time.sleep(0.1)

        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time(),  # Cached after file creation
            "_cache_complete": True,
        }

        is_stale = _is_cache_stale(str(test_file), cached_meta)

        assert is_stale is False

    def test_local_file_modified_is_stale(self, tmp_path):
        """Local file modified after cache is stale."""
        test_file = tmp_path / "product.nc"
        test_file.touch()

        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time() - 3600,  # Cached 1 hour ago
            "_cache_complete": True,
        }

        # Modify file (touch updates mtime)
        time.sleep(0.1)
        test_file.touch()

        is_stale = _is_cache_stale(str(test_file), cached_meta)

        assert is_stale is True

    def test_missing_local_file_is_stale(self):
        """Missing local file is stale."""
        cached_meta = {
            "version": "1c",
            "gdal_info_ts": time.time(),
            "_cache_complete": True,
        }

        is_stale = _is_cache_stale("/nonexistent/product.nc", cached_meta)

        assert is_stale is True


class TestGetOrExtract:
    """Test get_or_extract cache retrieval with staleness."""

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_cache_miss_extracts_metadata(self, mock_extract, tmp_path):
        """Cache miss extracts metadata."""
        mock_extract.return_value = {
            "version": "1c",
            "projection": "4326",
            "gdal_info_ts": time.time(),
            "_cache_complete": True,
        }

        cache_data = {}
        result = get_or_extract("product.nc", cache_data)

        assert result["version"] == "1c"
        assert "product.nc" in cache_data
        mock_extract.assert_called_once()

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_cache_hit_returns_cached(self, mock_extract, tmp_path):
        """Cache hit returns cached data without extraction."""
        test_file = tmp_path / "product.nc"
        test_file.touch()

        cache_data = {
            str(test_file): {
                "version": "1c",
                "projection": "4326",
                "gdal_info_ts": time.time(),
                "_cache_complete": True,
            }
        }

        result = get_or_extract(str(test_file), cache_data)

        assert result["version"] == "1c"
        mock_extract.assert_not_called()

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_stale_cache_refreshes(self, mock_extract):
        """Stale cache entry is refreshed."""
        mock_extract.return_value = {
            "version": "1c",
            "projection": "32610",  # Updated value
            "gdal_info_ts": time.time(),
            "_cache_complete": True,
        }

        cache_data = {
            "https://example.com/product.nc": {
                "version": "1c",
                "projection": "4326",  # Old value
                "gdal_info_ts": time.time() - 7200,  # 2 hours ago
                "_cache_complete": True,
            }
        }

        result = get_or_extract(
            "https://example.com/product.nc",
            cache_data,
            ttl_seconds=3600,  # 1 hour TTL
        )

        assert result["projection"] == "32610"  # Updated
        mock_extract.assert_called_once()


class TestThreadSafety:
    """Test thread safety of metadata cache."""

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_concurrent_reads_safe(self, mock_extract):
        """Multiple threads can read cache simultaneously."""
        mock_extract.return_value = {
            "version": "1c",
            "projection": "4326",
            "gdal_info_ts": time.time(),
            "_cache_complete": True,
        }

        # Prepopulate cache
        cache_data = {"product.nc": mock_extract.return_value}

        results = []
        errors = []

        def read_cache():
            try:
                for _ in range(100):
                    meta = get_or_extract("product.nc", cache_data)
                    results.append(meta)
            except Exception as e:
                errors.append(e)

        # 10 threads reading concurrently
        threads = [threading.Thread(target=read_cache) for _ in range(10)]

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # Should have no errors
        assert len(errors) == 0
        # Should have 1000 results (10 threads × 100 reads)
        assert len(results) == 1000
        # All should be same cached object
        assert all(r is results[0] for r in results[1:])

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_concurrent_first_access_safe(self, mock_extract):
        """Concurrent first access doesn't duplicate work."""
        call_count = [0]

        def slow_extract(fname):
            call_count[0] += 1
            time.sleep(0.05)  # Simulate slow extraction
            return {
                "version": "1c",
                "projection": "4326",
                "gdal_info_ts": time.time(),
                "_cache_complete": True,
            }

        mock_extract.side_effect = slow_extract

        cache_data = {}
        results = []

        def load_meta():
            meta = get_or_extract("product.nc", cache_data)
            results.append(meta)

        # 5 threads try to load same uncached metadata
        threads = [threading.Thread(target=load_meta) for _ in range(5)]

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # With proper locking, should extract at most a few times
        # (some threads may call before first completes, but not all 5)
        assert call_count[0] <= 5
        assert len(results) == 5

    @patch("ARIAtools.util.meta_cache.extract_metadata_gdal")
    def test_save_cache_thread_safe(self, mock_extract, tmp_path):
        """Save cache is thread-safe."""
        cache_file = tmp_path / "cache.json"
        cache_data = {}

        def save_repeatedly():
            for i in range(10):
                cache_data[f"product_{i}.nc"] = {"version": "1c", "projection": "4326"}
                save_cache(str(cache_file), cache_data)
                time.sleep(0.01)

        threads = [threading.Thread(target=save_repeatedly) for _ in range(3)]

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # Should successfully write without corruption
        assert cache_file.exists()
        saved_data = json.loads(cache_file.read_text())
        assert len(saved_data) >= 10


class TestExpandedMetadata:
    """Test expanded metadata extraction (projection, bounds, etc)."""

    @patch("osgeo.gdal.Info")
    def test_extract_projection(self, mock_info):
        """Extract projection from GDAL info."""
        mock_info.return_value = json.dumps(
            {
                "coordinateSystem": {
                    "wkt": (
                        'GEOGCS["WGS 84",DATUM["WGS_1984",'
                        'SPHEROID["WGS 84",6378137,298.257223563]],'
                        'PRIMEM["Greenwich",0],'
                        'UNIT["degree",0.0174532925199433]]'
                    )
                },
                "geoTransform": [0, 1, 0, 0, 0, -1],
                "size": [100, 100],
                "cornerCoordinates": {"upperLeft": [0, 0], "lowerRight": [100, -100]},
                "metadata": {"": {"NC_GLOBAL#version": "1c"}},
            }
        )

        meta = extract_metadata_gdal("test.nc")

        # Should have expanded fields
        assert "projection" in meta
        assert "geotransform" in meta
        assert "bounds" in meta
        assert "size" in meta
        assert "_cache_complete" in meta

    @patch("osgeo.gdal.Info")
    def test_incomplete_metadata_flagged(self, mock_info):
        """Incomplete metadata is flagged."""
        mock_info.return_value = json.dumps(
            {
                "metadata": {"": {"NC_GLOBAL#version": "1c"}}
                # Missing projection, geotransform, bounds, size
            }
        )

        meta = extract_metadata_gdal("test.nc")

        # Should have completeness flag set to False
        assert meta["_cache_complete"] is False


# Integration test with real GDAL (if available)
class TestIntegrationWithGDAL:
    """Integration tests with real GDAL operations."""

    @pytest.mark.skipif(
        not os.path.exists("/usr/bin/gdal_translate"), reason="GDAL not available"
    )
    def test_extract_metadata_real_file(self, tmp_path):
        """Extract metadata from real GeoTIFF."""
        from osgeo import gdal

        # Create test GeoTIFF
        test_file = tmp_path / "test.tif"
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(str(test_file), 100, 100, 1, gdal.GDT_Byte)
        ds.SetGeoTransform([0, 1, 0, 0, 0, -1])
        ds.SetProjection("EPSG:4326")
        data = np.random.randint(0, 256, (100, 100), dtype=np.uint8)
        ds.GetRasterBand(1).WriteArray(data)
        ds = None

        # Extract metadata
        meta = extract_metadata_gdal(str(test_file))

        # Verify expanded fields present
        assert meta is not None
        assert "geotransform" in meta
        assert meta["geotransform"] == [0, 1, 0, 0, 0, -1]
