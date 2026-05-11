"""Unit tests for the modern extract/timeseries command wrappers."""

from __future__ import annotations

import builtins
import sys
from types import ModuleType

import pytest

from aria_tools.commands import extract, timeseries


@pytest.mark.parametrize(
    ("module", "legacy_module_name", "expected_argv"),
    [
        (
            extract,
            "aria_tools.commands.legacy.extract",
            ["aria-tools extract", "--flag", "value"],
        ),
        (
            timeseries,
            "aria_tools.commands.legacy.timeseries",
            ["aria-tools timeseries", "--flag", "value"],
        ),
    ],
)
def test_wrapper_forwards_argv_and_restores_sys_argv(
    monkeypatch: pytest.MonkeyPatch,
    module: ModuleType,
    legacy_module_name: str,
    expected_argv: list[str],
) -> None:
    original_argv = ["python", "existing"]
    captured: list[str] | None = None

    legacy_module = ModuleType("legacy")

    def fake_main() -> None:
        nonlocal captured
        captured = sys.argv[:]

    legacy_module.main = fake_main  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, legacy_module_name, legacy_module)
    monkeypatch.setattr(sys, "argv", original_argv[:])

    module.main(["--flag", "value"])

    assert captured == expected_argv
    assert sys.argv == original_argv


@pytest.mark.parametrize(
    ("module", "missing_name", "script_name"),
    [
        (extract, "tile_mate", "ariaExtract.py"),
        (timeseries, "osgeo", "ariaTSsetup.py"),
    ],
)
def test_wrapper_prints_static_help_when_optional_dependency_is_missing(
    monkeypatch: pytest.MonkeyPatch,
    module: ModuleType,
    missing_name: str,
    script_name: str,
) -> None:
    printed: list[str] = []
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.startswith("aria_tools.commands.legacy."):
            raise ModuleNotFoundError(
                f"No module named '{missing_name}'",
                name=missing_name,
            )
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    monkeypatch.setattr(module, "_print_static_help", printed.append)

    module.main(["--help"])

    assert printed == [script_name]


@pytest.mark.parametrize("module", [extract, timeseries])
def test_wrapper_does_not_mask_real_import_errors(
    monkeypatch: pytest.MonkeyPatch,
    module: ModuleType,
) -> None:
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.startswith("aria_tools.commands.legacy."):
            raise ModuleNotFoundError(
                "No module named 'aria_tools.missing'",
                name="aria_tools.missing",
            )
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ModuleNotFoundError):
        module.main(["--help"])
