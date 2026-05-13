"""Week 4 tests for ``ariaOrderASF`` missing-dependency guidance."""

from __future__ import annotations

import builtins
import datetime as dt
import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "bin" / "ariaOrderASF.py"


def load_order_module(
    monkeypatch: pytest.MonkeyPatch,
    *,
    module_name: str,
    block_matplotlib: bool = False,
):
    """Load ``ariaOrderASF.py`` from disk with optional matplotlib blocking."""
    original_import = builtins.__import__

    if block_matplotlib:

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "matplotlib" or name.startswith("matplotlib."):
                raise ImportError("No module named 'matplotlib'")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", guarded_import)

    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_order_module_can_build_parser_without_matplotlib(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_order_module(
        monkeypatch,
        module_name="aria_order_no_matplotlib_parser",
        block_matplotlib=True,
    )

    parser = module.create_parser()
    args = parser.parse_args(["--statusjobs", "--status-name", "example-job"])

    assert args.statusjobs is True
    assert args.status_name == "example-job"


def test_order_plotting_points_to_install_when_matplotlib_is_missing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_order_module(
        monkeypatch,
        module_name="aria_order_no_matplotlib_plotting",
        block_matplotlib=True,
    )

    with pytest.raises(ImportError, match=r"pip install -e \."):
        module.plot_baseline(
            {
                dt.date(2024, 1, 1): 0.0,
                dt.date(2024, 2, 1): 1.5,
            },
            {"existing": [], "new": []},
            25050,
            output_dir=str(tmp_path),
        )


def test_order_missing_dependency_errors_reference_full_install(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_order_module(
        monkeypatch,
        module_name="aria_order_optional_dependency_messages",
    )

    monkeypatch.setattr(module, "HAS_ASF_ENUMERATION", False)
    with pytest.raises(ImportError, match=r"pip install -e \."):
        module.get_acquisitions_for_frame(
            25050,
            dt.date(2024, 1, 1),
            dt.date(2024, 2, 1),
        )

    monkeypatch.setattr(module, "HAS_HYP3_SDK", False)
    with pytest.raises(ImportError, match=r"pip install -e \."):
        module.order_pairs(25050, "unused.csv")
