"""Unit tests for shared CLI helpers."""

from __future__ import annotations

from aria_tools.cli.common import inject_workdir


def test_inject_workdir_adds_missing_workdir() -> None:
    assert inject_workdir(["--flag", "value"], "/tmp/job") == [
        "--workdir",
        "/tmp/job",
        "--flag",
        "value",
    ]


def test_inject_workdir_preserves_existing_long_flag() -> None:
    assert inject_workdir(["--workdir", "/tmp/existing"], "/tmp/job") == [
        "--workdir",
        "/tmp/existing",
    ]


def test_inject_workdir_preserves_existing_short_flag() -> None:
    assert inject_workdir(["-w", "/tmp/existing"], "/tmp/job") == [
        "-w",
        "/tmp/existing",
    ]


def test_inject_workdir_is_noop_without_top_level_value() -> None:
    assert inject_workdir(["--flag"], None) == ["--flag"]
