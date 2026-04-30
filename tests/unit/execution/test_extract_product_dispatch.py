"""Unit tests for export backend selection in ``ARIAtools.extractProduct``."""

from __future__ import annotations

import ARIAtools.extractProduct as extract_product


class FakeProgressBar:
    def __init__(self, *, maxValue: int, prefix: str) -> None:
        self.max_value = maxValue
        self.prefix = prefix
        self.updates: list[int] = []
        self.closed = False

    def update(self, completed: int) -> None:
        self.updates.append(completed)

    def close(self) -> None:
        self.closed = True


def test_run_export_jobs_uses_serial_backend_for_single_mode(
    monkeypatch,
) -> None:
    progress_bar = FakeProgressBar(maxValue=2, prefix="Exporting layer: ")

    monkeypatch.setattr(
        extract_product.ARIAtools.util.misc,
        "ProgressBar",
        lambda **kwargs: progress_bar,
    )
    monkeypatch.setattr(
        extract_product,
        "export_product_worker_helper",
        lambda arg: f"serial-{arg}",
    )

    outputs = extract_product._run_export_jobs(
        ["a", "b"], num_workers=4, layer="layer", multiproc_method="single"
    )

    assert outputs == ["serial-a", "serial-b"]
    assert progress_bar.updates == [1, 2]
    assert progress_bar.closed is True


def test_run_export_jobs_falls_back_to_serial_for_single_worker(
    monkeypatch,
) -> None:
    progress_bar = FakeProgressBar(maxValue=1, prefix="Exporting layer: ")

    monkeypatch.setattr(
        extract_product.ARIAtools.util.misc,
        "ProgressBar",
        lambda **kwargs: progress_bar,
    )
    monkeypatch.setattr(
        extract_product,
        "export_product_worker_helper",
        lambda arg: f"serial-{arg}",
    )
    monkeypatch.setattr(
        extract_product,
        "_run_export_products_with_processes",
        lambda *args, **kwargs: ["processes-should-not-run"],
    )

    outputs = extract_product._run_export_jobs(
        ["only"], num_workers=1, layer="layer", multiproc_method="processes"
    )

    assert outputs == ["serial-only"]
    assert progress_bar.updates == [1]
    assert progress_bar.closed is True


def test_run_export_jobs_uses_process_backend_when_available(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_run_processes(mp_args, num_workers, layer):
        captured["mp_args"] = mp_args
        captured["num_workers"] = num_workers
        captured["layer"] = layer
        return ["process-output"]

    monkeypatch.setattr(
        extract_product,
        "_run_export_products_with_processes",
        fake_run_processes,
    )

    outputs = extract_product._run_export_jobs(
        ["a"], num_workers=3, layer="layer", multiproc_method="processes"
    )

    assert outputs == ["process-output"]
    assert captured == {"mp_args": ["a"], "num_workers": 3, "layer": "layer"}


def test_run_export_jobs_rejects_removed_thread_backend() -> None:
    try:
        extract_product._run_export_jobs(
            ["a"], num_workers=2, layer="layer", multiproc_method="threads"
        )
    except ValueError as exc:
        assert '"single" or "processes"' in str(exc)
    else:
        raise AssertionError("Expected ValueError for removed thread backend")
