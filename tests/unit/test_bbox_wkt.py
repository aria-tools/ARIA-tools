"""Unit tests for shared bbox parsing, including WKT polygons."""

from __future__ import annotations

from pathlib import Path

import pytest

from ARIAtools import product


def test_parse_bbox_argument_accepts_wkt_polygon(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def fake_save_shp(fname, polygon, projection, drivername="GeoJSON"):
        captured["fname"] = fname
        captured["polygon"] = polygon
        captured["projection"] = projection
        captured["drivername"] = drivername

    monkeypatch.setattr(product.ARIAtools.util.shp, "save_shp", fake_save_shp)

    polygon, bbox_file = product._parse_bbox_argument(
        "POLYGON((-118 36, -118 37, -117 37, -117 36, -118 36))",
        str(tmp_path),
        4326,
    )

    assert polygon.bounds == (-118.0, 36.0, -117.0, 37.0)
    assert bbox_file == str(tmp_path / "user_bbox.json")
    assert captured == {
        "fname": str(tmp_path / "user_bbox.json"),
        "polygon": polygon,
        "projection": 4326,
        "drivername": "GeoJSON",
    }


def test_parse_bbox_argument_accepts_snwe_string(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        product.ARIAtools.util.shp,
        "save_shp",
        lambda *args, **kwargs: None,
    )

    polygon, bbox_file = product._parse_bbox_argument(
        "36 37 -118 -117",
        str(tmp_path),
        4326,
    )

    assert polygon.bounds == (-118.0, 36.0, -117.0, 37.0)
    assert bbox_file == str(tmp_path / "user_bbox.json")


def test_parse_bbox_argument_rejects_invalid_text(tmp_path: Path) -> None:
    with pytest.raises(Exception, match="Cannot understand the --bbox argument"):
        product._parse_bbox_argument("not-a-bbox", str(tmp_path), 4326)

