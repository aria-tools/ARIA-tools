"""Top-level ``aria-tools`` command router."""

from __future__ import annotations

import argparse
import importlib
import sys
from collections.abc import Callable, Sequence

from aria_tools import __version__
from aria_tools.cli.common import inject_workdir
from aria_tools.config.logging import (
    DEFAULT_LOG_LEVEL,
    LOG_LEVELS,
    configure_logging,
)
from aria_tools.errors import AriaToolsError

COMMANDS = {
    "download": "aria_tools.commands.download",
    "extract": "aria_tools.commands.extract",
    "timeseries": "aria_tools.commands.timeseries",
    "plot": "aria_tools.commands.plot",
    "order": "aria_tools.commands.order",
    "misclosure": "aria_tools.commands.misclosure",
    "aoi": "aria_tools.commands.aoi",
    "kml2box": "aria_tools.commands.kml",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aria-tools",
        description="Modern command shell for ARIA-tools workflows.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--log-level",
        choices=LOG_LEVELS,
        default=DEFAULT_LOG_LEVEL,
        help="Logger log level for the command shell. Default: info.",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        default=None,
        help=(
            "Shared working directory. Forwarded to subcommands unless the "
            "subcommand arguments already include -w/--workdir."
        ),
    )

    command_list = ", ".join(COMMANDS)
    parser.add_argument("command", choices=COMMANDS, help=f"One of: {command_list}")
    parser.add_argument(
        "command_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed through to the selected command.",
    )
    return parser


def _load_command(command: str) -> Callable[[Sequence[str] | None], None]:
    module = importlib.import_module(COMMANDS[command])
    return module.main


def run(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)

    command_args = inject_workdir(args.command_args, args.workdir)
    command_main = _load_command(args.command)
    command_main(command_args)
    return 0


def main(argv: Sequence[str] | None = None) -> None:
    try:
        raise SystemExit(run(argv))
    except AriaToolsError as exc:
        print(f"aria-tools: error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
