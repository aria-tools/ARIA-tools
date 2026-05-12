"""Unit tests for metadata cache helpers."""

from __future__ import annotations

import json
from pathlib import Path

from ARIAtools.util import meta_cache


def test_cache_path_uses_url_file_directory(tmp_path: Path) -> None:
    url_file = tmp_path / "products.txt"
    url_file.write_text("https://example.test/file.nc\n")

    assert meta_cache._cache_path(str(url_file)) == str(
        tmp_path / "aria_meta_cache.json"
    )


def test_file_key_strips_virtual_filesystem_prefixes() -> None:
    assert (
        meta_cache._file_key("/vsicurl/https://example.test/file.nc")
        == "https://example.test/file.nc"
    )
    assert meta_cache._file_key("/vsis3/example-bucket/path/file.nc") == (
        "example-bucket/path/file.nc"
    )


def test_save_cache_and_load_cache_round_trip(tmp_path: Path) -> None:
    cache_file = tmp_path / "aria_meta_cache.json"
    cache_data = {"https://example.test/file.nc": {"version": "2_0_0"}}

    meta_cache.save_cache(str(cache_file), cache_data)

    assert json.loads(cache_file.read_text()) == cache_data
    assert meta_cache.load_cache(str(cache_file)) == cache_data


def test_load_cache_returns_empty_dict_for_invalid_json(tmp_path: Path) -> None:
    cache_file = tmp_path / "aria_meta_cache.json"
    cache_file.write_text("{not-json")

    assert meta_cache.load_cache(str(cache_file)) == {}


def test_get_or_extract_uses_cached_entry_without_reextracting(
    monkeypatch,
) -> None:
    cache_data = {"https://example.test/file.nc": {"version": "cached"}}
    calls: list[str] = []
    monkeypatch.setattr(
        meta_cache,
        "extract_metadata_gdal",
        lambda fname: calls.append(fname),
    )

    result = meta_cache.get_or_extract(
        "/vsicurl/https://example.test/file.nc",
        cache_data,
    )

    assert result == {"version": "cached"}
    assert calls == []


def test_get_or_extract_populates_cache_on_miss(monkeypatch) -> None:
    extracted = {"version": "2_0_1", "driver": "netCDF"}
    cache_data: dict[str, dict[str, str]] = {}
    monkeypatch.setattr(
        meta_cache,
        "extract_metadata_gdal",
        lambda fname: extracted,
    )

    result = meta_cache.get_or_extract(
        "/vsicurl/https://example.test/file.nc",
        cache_data,
    )

    assert result == extracted
    assert cache_data == {"https://example.test/file.nc": extracted}
