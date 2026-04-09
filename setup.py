#!/usr/bin/env python3
"""Compatibility shim for tools that still invoke ``python setup.py``.

Project metadata and install behavior live in ``pyproject.toml``.
"""

from setuptools import setup


setup()
