"""Integration tests for clear CLI failures on missing input products."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
TOOLS_ROOT = REPO_ROOT / "tools"


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    pythonpath_parts = [str(SRC_ROOT), str(TOOLS_ROOT)]
    if env.get("PYTHONPATH"):
        pythonpath_parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    return subprocess.run(
        [sys.executable, "-m", "aria_tools.cli.app", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env=env,
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
        "No input products matched the provided file argument: "
        "does_not_exist/*.nc"
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
        "No input products matched the provided file argument: "
        "does_not_exist/*.nc"
    ) in result.stderr
    assert "Traceback" not in result.stderr
