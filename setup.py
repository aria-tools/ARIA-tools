#!/usr/bin/env python3
"""Compatibility shim for tools that still invoke ``python setup.py``.

Project metadata and install behavior live in ``pyproject.toml``.
"""


def main() -> None:
    from setuptools import setup

    setup()


if __name__ == "__main__":
    main()
