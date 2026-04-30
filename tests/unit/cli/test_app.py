"""Unit tests for the modern ``aria-tools`` router."""

from __future__ import annotations

import types

import pytest

from aria_tools import __version__
from aria_tools.cli import app
from aria_tools.config.logging import DEFAULT_LOG_LEVEL
from aria_tools.errors import AriaToolsError


def test_build_parser_exposes_expected_defaults() -> None:
    parser = app.build_parser()

    args = parser.parse_args(["download"])
    command_action = next(
        action for action in parser._actions if action.dest == "command"
    )

    assert parser.prog == "aria-tools"
    assert args.log_level == DEFAULT_LOG_LEVEL
    assert args.workdir is None
    assert args.command == "download"
    assert args.command_args == []
    assert tuple(command_action.choices) == tuple(app.COMMANDS)


def test_main_prints_top_level_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        app.main(["--help"])

    captured = capsys.readouterr()

    assert excinfo.value.code == 0
    assert "Modern command shell for ARIA-tools workflows." in captured.out
    assert "download" in captured.out
    assert captured.err == ""


def test_main_prints_top_level_version(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        app.main(["--version"])

    captured = capsys.readouterr()

    assert excinfo.value.code == 0
    assert captured.out.strip() == f"aria-tools {__version__}"
    assert captured.err == ""


def test_run_routes_command_and_injects_workdir(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: dict[str, object] = {}

    def fake_main(argv: list[str] | None = None) -> None:
        calls["argv"] = argv

    def fake_import_module(module_name: str) -> types.SimpleNamespace:
        calls["module_name"] = module_name
        return types.SimpleNamespace(main=fake_main)

    log_levels: list[str] = []
    monkeypatch.setattr(app.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(app, "configure_logging", log_levels.append)

    exit_code = app.run(
        [
            "--log-level",
            "debug",
            "--workdir",
            "/tmp/job",
            "extract",
            "--flag",
            "value",
        ]
    )

    assert exit_code == 0
    assert log_levels == ["debug"]
    assert calls["module_name"] == "aria_tools.commands.extract"
    assert calls["argv"] == ["--workdir", "/tmp/job", "--flag", "value"]


def test_run_does_not_duplicate_subcommand_workdir(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_argv: list[str] | None = None

    def fake_main(argv: list[str] | None = None) -> None:
        nonlocal captured_argv
        captured_argv = argv

    monkeypatch.setattr(
        app.importlib,
        "import_module",
        lambda module_name: types.SimpleNamespace(main=fake_main),
    )
    monkeypatch.setattr(app, "configure_logging", lambda level: None)

    app.run(
        [
            "--workdir",
            "/tmp/top-level",
            "extract",
            "--workdir",
            "/tmp/subcommand",
        ]
    )

    assert captured_argv == ["--workdir", "/tmp/subcommand"]


def test_main_converts_expected_errors_to_exit_code_2(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def raise_expected_error(argv: object = None) -> int:
        raise AriaToolsError("bad inputs")

    monkeypatch.setattr(app, "run", raise_expected_error)

    with pytest.raises(SystemExit) as excinfo:
        app.main(["extract"])

    captured = capsys.readouterr()

    assert excinfo.value.code == 2
    assert captured.out == ""
    assert captured.err.strip() == "aria-tools: error: bad inputs"
