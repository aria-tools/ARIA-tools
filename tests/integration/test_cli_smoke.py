"""Smoke tests for the modern CLI entry point."""

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


def test_cli_help_smoke() -> None:
    result = run_cli("--help")

    assert result.returncode == 0
    assert "Modern command shell for ARIA-tools workflows." in result.stdout
    assert "download" in result.stdout
    assert result.stderr == ""


def test_cli_version_smoke() -> None:
    result = run_cli("--version")

    assert result.returncode == 0
    assert result.stdout.strip().startswith("aria-tools ")
    assert result.stderr == ""
