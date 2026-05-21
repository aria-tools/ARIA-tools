#! /usr/bin/env python3
"""Compatibility wrapper for the migrated extract command."""


def run() -> None:
    from aria_tools.commands.extract import main

    main()


if __name__ == "__main__":
    run()
