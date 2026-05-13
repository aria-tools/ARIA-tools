#! /usr/bin/env python3
"""Compatibility wrapper for the migrated time-series setup command."""

import sys
from pathlib import Path

src_dir = Path(__file__).resolve().parents[2] / "src"
if src_dir.exists():
    sys.path.insert(0, str(src_dir))


def run() -> None:
    from aria_tools.commands.timeseries import main

    main()


if __name__ == "__main__":
    run()
