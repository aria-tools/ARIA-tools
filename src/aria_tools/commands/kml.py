"""KML-to-bounding-box command adapter."""

from __future__ import annotations

from aria_tools.cli.types import Argv
from aria_tools.commands._legacy import run_legacy_script


def main(argv: Argv = None) -> None:
    run_legacy_script("ariaKml2box.py", argv)
