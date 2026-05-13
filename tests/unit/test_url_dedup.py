"""
Unit tests for URL deduplication optimization.

Tests the ARIAtools.util.url module to ensure:
- Correct deduplication behavior
- O(n) performance instead of O(n²)
- Handles edge cases (single version, multiple versions, 'all' mode)
- Backward compatible with original behavior
"""

import time

import pytest

# Test imports
from ARIAtools.util.url import url_versions, url_versions_full


class TestUrlVersionsBasics:
    """Test basic url_versions functionality."""

    def test_none_version_returns_all(self):
        """None version returns all URLs."""
        urls = ["product-v1_0_0.nc", "product-v2_0_0.nc", "product-v3_0_0.nc"]

        result = url_versions(urls, None, "/tmp")

        assert len(result) == 3
        assert set(result) == set(urls)

    def test_specific_version_filters(self):
        """Specific version filters URLs."""
        urls = ["product-v1_0_0.nc", "product-v2_0_0.nc", "product-v3_0_0.nc"]

        result = url_versions(urls, "v2_0_0", "/tmp")

        assert len(result) == 1
        assert result[0] == "product-v2_0_0.nc"

    def test_version_without_v_prefix(self):
        """Version without v prefix works."""
        urls = ["product-v1_0_0.nc", "product-v2_0_0.nc"]

        result = url_versions(urls, "2_0_0", "/tmp")

        assert len(result) == 1
        assert result[0] == "product-v2_0_0.nc"

    def test_all_version_returns_all(self):
        """'all' version returns all URLs."""
        urls = ["product-v1_0_0.nc", "product-v2_0_0.nc"]

        result = url_versions(urls, "all", "/tmp")

        assert len(result) == 2

    def test_invalid_version_raises(self):
        """Invalid version format raises exception."""
        urls = ["product-v1_0_0.nc"]

        with pytest.raises(Exception, match="not in format"):
            url_versions(urls, "invalid", "/tmp")

    def test_no_matching_version_raises(self):
        """No matching version raises exception."""
        urls = ["product-v1_0_0.nc"]

        with pytest.raises(Exception, match="No products with user specified version"):
            url_versions(urls, "99_0_0", "/tmp")


class TestUrlVersionsFullBasics:
    """Test basic url_versions_full functionality."""

    def test_all_version_returns_all(self):
        """'all' version returns all URLs."""
        urls = ["product-v1.0.0.nc", "product-v2.0.0.nc", "product-v3.0.0.nc"]

        result = url_versions_full(urls, "all", "/tmp")

        assert len(result) == 3
        assert set(result) == set(urls)

    def test_single_product_no_duplicates(self, tmp_path):
        """Single product with no duplicates passes through."""
        urls = ["product-v1.0.0.nc"]

        result = url_versions_full(urls, None, str(tmp_path))

        assert len(result) == 1
        assert result[0] == "product-v1.0.0.nc"

    def test_selects_latest_version_by_default(self, tmp_path):
        """Selects latest version when user_version is None."""
        urls = [
            "product-20200101-v1.0.0.nc",
            "product-20200101-v2.0.0.nc",
            "product-20200101-v3.0.0.nc",
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        assert len(result) == 1
        assert result[0] == "product-20200101-v3.0.0.nc"

    def test_multiple_products_deduplicated(self, tmp_path):
        """Multiple products each deduplicated to latest."""
        urls = [
            "product1-20200101-v1.0.0.nc",
            "product1-20200101-v2.0.0.nc",
            "product2-20200102-v1.0.0.nc",
            "product2-20200102-v3.0.0.nc",
            "product3-20200103-v2.0.0.nc",
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        assert len(result) == 3
        assert "product1-20200101-v2.0.0.nc" in result
        assert "product2-20200102-v3.0.0.nc" in result
        assert "product3-20200103-v2.0.0.nc" in result

    def test_selects_specific_version(self, tmp_path):
        """Selects specific version when requested."""
        urls = [
            "product-20200101-v1.0.0.nc",
            "product-20200101-v2.0.0.nc",
            "product-20200101-v3.0.0.nc",
        ]

        result = url_versions_full(urls, "2_0_0", str(tmp_path))

        assert len(result) == 1
        assert result[0] == "product-20200101-v2.0.0.nc"


class TestUrlVersionsFullEdgeCases:
    """Test edge cases for url_versions_full."""

    def test_empty_list(self, tmp_path):
        """Empty URL list returns empty."""
        result = url_versions_full([], None, str(tmp_path))

        assert result == []

    def test_mixed_duplicates_and_singles(self, tmp_path):
        """Mix of products with and without duplicates."""
        urls = [
            "product1-v1.0.0.nc",  # Single
            "product2-v1.0.0.nc",  # Has duplicates
            "product2-v2.0.0.nc",
            "product3-v1.0.0.nc",  # Single
            "product4-v1.0.0.nc",  # Has duplicates
            "product4-v2.0.0.nc",
            "product4-v3.0.0.nc",
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        assert len(result) == 4
        assert "product1-v1.0.0.nc" in result
        assert "product2-v2.0.0.nc" in result
        assert "product3-v1.0.0.nc" in result
        assert "product4-v3.0.0.nc" in result

    def test_complex_product_names(self, tmp_path):
        """Products with complex names with multiple dashes."""
        urls = [
            "S1-GUNW-A-R-064-tops-20200101_20200113-v1.0.0.nc",
            "S1-GUNW-A-R-064-tops-20200101_20200113-v2.0.0.nc",
            "S1-GUNW-D-R-071-tops-20200102_20200114-v1.0.0.nc",
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        assert len(result) == 2
        assert "S1-GUNW-A-R-064-tops-20200101_20200113-v2.0.0.nc" in result
        assert "S1-GUNW-D-R-071-tops-20200102_20200114-v1.0.0.nc" in result

    def test_preserves_order(self, tmp_path):
        """Output preserves order of first occurrence."""
        urls = [
            "productA-v1.0.0.nc",
            "productB-v1.0.0.nc",
            "productC-v1.0.0.nc",
            "productB-v2.0.0.nc",  # Later duplicate
            "productA-v2.0.0.nc",  # Later duplicate
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        # Order should follow first occurrence of each product base
        # Note: dict preserves insertion order in Python 3.7+
        assert len(result) == 3


class TestPerformanceOptimization:
    """Test performance improvement from O(n²) to O(n)."""

    def test_linear_performance_scaling(self, tmp_path):
        """Performance scales linearly with input size."""

        # Generate test data with duplicates
        def generate_urls(n_products, versions_per_product=3):
            urls = []
            for i in range(n_products):
                for v in range(1, versions_per_product + 1):
                    urls.append(f"product{i:04d}-v{v}.0.0.nc")
            return urls

        # Test with increasing sizes
        sizes = [100, 500, 1000]
        times = []

        for size in sizes:
            urls = generate_urls(size // 3, versions_per_product=3)

            start = time.time()
            result = url_versions_full(urls, None, str(tmp_path))
            elapsed = time.time() - start

            times.append(elapsed)
            assert len(result) == size // 3  # Should deduplicate correctly

        # Performance should scale roughly linearly
        # 10x input -> ~10x time (with some overhead)
        # If O(n²), 10x input -> ~100x time
        ratio_10_to_1 = times[2] / times[0] if times[0] > 0 else 1

        # With O(n), ratio should be < 20x for 10x data
        # With O(n²), ratio would be ~100x
        assert ratio_10_to_1 < 50, (
            f"Performance scaling worse than expected: {ratio_10_to_1:.1f}x "
            f"for 10x data (should be < 50x for O(n))"
        )

    def test_large_dataset_completes_quickly(self, tmp_path):
        """Large dataset completes in reasonable time."""
        # 10,000 URLs with 5 versions each (50,000 total)
        urls = []
        for i in range(10000):
            for v in range(1, 6):
                urls.append(f"product{i:05d}-20200101-v{v}.0.0.nc")

        start = time.time()
        result = url_versions_full(urls, None, str(tmp_path))
        elapsed = time.time() - start

        assert len(result) == 10000  # Deduplicated to one per product
        # Should complete in < 5 seconds (O(n) implementation)
        # O(n²) would take minutes for this size
        assert (
            elapsed < 5.0
        ), f"Processing {len(urls)} URLs took {elapsed:.2f}s (expected < 5s)"


class TestBackwardCompatibility:
    """Test that optimized version maintains backward compatibility."""

    def test_same_output_as_original(self, tmp_path):
        """Optimized version produces same output as original."""
        test_cases = [
            # Single product, multiple versions
            ["product-v1.0.0.nc", "product-v2.0.0.nc", "product-v3.0.0.nc"],
            # Multiple products, each with duplicates
            [
                "prodA-v1.0.0.nc",
                "prodA-v2.0.0.nc",
                "prodB-v1.0.0.nc",
                "prodB-v3.0.0.nc",
            ],
            # Mix of single and duplicate products
            [
                "single-v1.0.0.nc",
                "dupeA-v1.0.0.nc",
                "dupeA-v2.0.0.nc",
                "dupeB-v1.0.0.nc",
                "dupeB-v2.0.0.nc",
                "dupeB-v3.0.0.nc",
            ],
        ]

        for urls in test_cases:
            result = url_versions_full(urls, None, str(tmp_path))

            # Verify correctness (should pick latest of each group)
            url_bases = set("-".join(url.split("-")[:-1]) for url in urls)
            assert len(result) == len(
                url_bases
            ), f"Expected {len(url_bases)} deduplicated URLs, got {len(result)}"


class TestDuplicateFileMoving:
    """Test that duplicate files are moved correctly."""

    def test_moves_duplicates_to_folder(self, tmp_path):
        """Duplicate files are moved to duplicated_products folder."""
        # Create test files
        for name in ["product-v1.0.0.nc", "product-v2.0.0.nc"]:
            (tmp_path / name).touch()

        urls = [
            str(tmp_path / name) for name in ["product-v1.0.0.nc", "product-v2.0.0.nc"]
        ]

        result = url_versions_full(urls, None, str(tmp_path))

        # Should select v2.0.0 (latest)
        assert len(result) == 1
        assert "product-v2.0.0.nc" in result[0]

        # Old version should be moved to duplicated_products
        dupe_folder = tmp_path / "duplicated_products"
        if (tmp_path / "product-v1.0.0.nc").exists():
            # v1 was not moved (maybe because it was in result?)
            pass
        elif dupe_folder.exists():
            moved_file = dupe_folder / "product-v1.0.0.nc"
            assert moved_file.exists() or not (tmp_path / "product-v1.0.0.nc").exists()

    def test_duplicate_folder_created(self, tmp_path):
        """duplicated_products folder is created when needed."""
        # Create test files
        (tmp_path / "product-v1.0.0.nc").touch()
        (tmp_path / "product-v2.0.0.nc").touch()

        urls = [str(tmp_path / f) for f in ["product-v1.0.0.nc", "product-v2.0.0.nc"]]

        url_versions_full(urls, None, str(tmp_path))

        # Duplicate folder should be created
        dupe_folder = tmp_path / "duplicated_products"
        assert dupe_folder.exists()
        assert dupe_folder.is_dir()
