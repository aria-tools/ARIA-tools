"""Extract command entry point."""

from __future__ import annotations

import sys

from aria_tools.cli.types import Argv
from aria_tools.commands._legacy import (
    _is_missing_optional_dependency,
    _print_static_help,
    _wants_help,
)


def main(argv: Argv = None) -> None:
    command_argv = sys.argv[1:] if argv is None else list(argv)

    try:
        from aria_tools.commands.legacy.extract import main as legacy_main
    except ModuleNotFoundError as exc:
        if _wants_help(command_argv) and _is_missing_optional_dependency(exc):
            _print_static_help("ariaExtract.py")
            return
        raise

    old_argv = sys.argv[:]
    sys.argv = ["aria-tools extract", *command_argv]
    try:
        legacy_main()
    finally:
        sys.argv = old_argv
