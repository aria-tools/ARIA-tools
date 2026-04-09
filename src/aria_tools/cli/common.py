"""Shared CLI helpers for the modern ARIA-tools command router."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from aria_tools.errors import CommandDispatchError

REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_BIN_DIR = REPO_ROOT / "tools" / "bin"


def legacy_script_path(script_name: str) -> Path:
    """Return the source-tree path for a legacy command script."""

    script_path = LEGACY_BIN_DIR / script_name
    if not script_path.exists():
        raise CommandDispatchError(
            f"Cannot locate legacy command script: {script_path}"
        )
    return script_path


def inject_workdir(argv: Sequence[str], workdir: str | None) -> list[str]:
    """Forward a top-level workdir to commands that did not receive one."""

    args = list(argv)
    if workdir is None:
        return args
    if "-w" in args or "--workdir" in args:
        return args
    return ["--workdir", workdir, *args]
