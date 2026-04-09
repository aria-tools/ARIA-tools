"""Logging helpers shared by modern ARIA-tools commands."""

from __future__ import annotations

import logging

LOG_LEVELS = ("debug", "info", "warning", "error")
DEFAULT_LOG_LEVEL = "info"


def configure_logging(level: str = DEFAULT_LOG_LEVEL) -> None:
    """Configure process-wide logging for the modern command shell."""

    numeric_level = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }[level]
    logging.basicConfig(level=numeric_level)
