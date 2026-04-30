"""Shared pytest configuration for the source-tree test layout."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

for import_root in (REPO_ROOT / "src", REPO_ROOT / "tools"):
    import_root_str = str(import_root)
    if import_root_str not in sys.path:
        sys.path.insert(0, import_root_str)
