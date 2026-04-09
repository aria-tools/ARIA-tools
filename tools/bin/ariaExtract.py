#! /usr/bin/env python3
"""Compatibility wrapper for the migrated extract command."""

import sys
from pathlib import Path

src_dir = Path(__file__).resolve().parents[2] / "src"
if src_dir.exists():
    sys.path.insert(0, str(src_dir))

from aria_tools.commands.extract import main


if __name__ == "__main__":
    main()
