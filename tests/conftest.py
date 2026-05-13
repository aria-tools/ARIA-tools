"""Shared pytest configuration for install-first test runs."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = _REPO_ROOT / "src"
_TOOLS_DIR = _REPO_ROOT / "tools"

for _path in (_SRC_DIR, _TOOLS_DIR):
    _path_str = str(_path)
    if _path_str not in sys.path:
        sys.path.insert(0, _path_str)


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run tests marked slow.",
    )
    parser.addoption(
        "--run-network",
        action="store_true",
        default=False,
        help="Run tests marked network_required.",
    )
    parser.addoption(
        "--run-credentialed",
        action="store_true",
        default=False,
        help="Run tests marked credentialed.",
    )


def pytest_collection_modifyitems(
    config: pytest.Config,
    items: list[pytest.Item],
) -> None:
    skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
    skip_network = pytest.mark.skip(reason="need --run-network option to run")
    skip_credentialed = pytest.mark.skip(reason="need --run-credentialed option to run")

    for item in items:
        item_path = Path(str(getattr(item, "path", item.fspath)))
        path_parts = item_path.parts

        if "tests" in path_parts and (
            "unit" in path_parts or "integration" in path_parts
        ):
            item.add_marker(pytest.mark.offline)

        if "regression" in path_parts or item_path.name == "test_virtual_access.py":
            item.add_marker(pytest.mark.slow)
            item.add_marker(pytest.mark.network_required)
            item.add_marker(pytest.mark.credentialed)

        if "slow" in item.keywords and not config.getoption("--run-slow"):
            item.add_marker(skip_slow)
        if "network_required" in item.keywords and not config.getoption(
            "--run-network"
        ):
            item.add_marker(skip_network)
        if "credentialed" in item.keywords and not config.getoption(
            "--run-credentialed"
        ):
            item.add_marker(skip_credentialed)
