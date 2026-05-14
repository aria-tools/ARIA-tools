"""Unit tests for extract command validation paths."""

from __future__ import annotations

import pytest

from aria_tools.cli import app


def test_extract_requires_dem_for_dem_dependent_layers(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = app.main(
        [
            "extract",
            "-f",
            "does_not_exist/*.nc",
            "-l",
            "bPerpendicular",
        ]
    )

    captured = capsys.readouterr()

    assert exit_code == 2
    assert "A valid DEM must be specified when extracting any of" in captured.err
    assert "bPerpendicular" in captured.err
    assert "Traceback" not in captured.err
