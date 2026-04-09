"""Shared exception types for the modern ARIA-tools CLI."""

from __future__ import annotations


class AriaToolsError(Exception):
    """Base class for expected ARIA-tools command failures."""


class CommandDispatchError(AriaToolsError):
    """Raised when a CLI command cannot be dispatched."""
