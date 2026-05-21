#! /usr/bin/env python3
"""Compatibility wrapper for the migrated time-series setup command."""


def run() -> None:
    from aria_tools.commands.timeseries import main

    main()


if __name__ == "__main__":
    run()
