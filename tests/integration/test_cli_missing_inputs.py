"""Integration tests for clear CLI failures on missing input products."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "aria_tools", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_extract_fails_cleanly_when_no_input_products_match() -> None:
    result = run_cli(
        "extract",
        "-f",
        "does_not_exist/*.nc",
        "-l",
        "coherence",
        "--log-level",
        "info",
    )

    assert result.returncode == 2
    assert (
        "No input products matched the provided file argument: does_not_exist/*.nc"
    ) in result.stderr
    assert "Traceback" not in result.stderr


def test_timeseries_fails_cleanly_when_no_input_products_match() -> None:
    result = run_cli(
        "timeseries",
        "-f",
        "does_not_exist/*.nc",
        "-l",
        "coherence",
        "--log-level",
        "info",
    )

    assert result.returncode == 2
    assert (
        "No input products matched the provided file argument: does_not_exist/*.nc"
    ) in result.stderr
    assert "Traceback" not in result.stderr
