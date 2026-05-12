"""Unit tests for native worker execution helpers."""

from __future__ import annotations

import pytest

from aria_tools.execution import workers
from aria_tools.execution.workers import resolve_worker_count, run_process_pool


def multiply_by_two(value: int) -> int:
    return value * 2


def fail_on_two(value: int) -> int:
    if value == 2:
        raise ValueError("worker boom")
    return value


class FakeFuture:
    def __init__(self, value: int, worker) -> None:
        self._value = value
        self._worker = worker

    def result(self) -> int:
        return self._worker(self._value)


class FakeProcessPoolExecutor:
    def __init__(self, *, max_workers: int) -> None:
        self.max_workers = max_workers
        self.submitted: list[FakeFuture] = []

    def __enter__(self) -> FakeProcessPoolExecutor:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def submit(self, worker, item: int) -> FakeFuture:
        future = FakeFuture(item, worker)
        self.submitted.append(future)
        return future


def test_resolve_worker_count_expands_all_cpus(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("aria_tools.execution.workers.os.cpu_count", lambda: 12)

    assert resolve_worker_count("ALL_CPUS") == 12


@pytest.mark.parametrize("raw_value, expected", [(3, 3), ("4", 4)])
def test_resolve_worker_count_coerces_numeric_values(
    raw_value: int | str,
    expected: int,
) -> None:
    assert resolve_worker_count(raw_value) == expected


def test_run_process_pool_collects_results_and_progress() -> None:
    completed: list[int] = []
    executor = FakeProcessPoolExecutor(max_workers=2)

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(
        workers.concurrent.futures,
        "ProcessPoolExecutor",
        lambda max_workers: executor,
    )
    monkeypatch.setattr(
        workers.concurrent.futures,
        "as_completed",
        lambda futures: reversed(list(futures)),
    )

    try:
        results = run_process_pool(
            multiply_by_two,
            [1, 2, 3],
            max_workers=2,
            on_complete=completed.append,
        )
    finally:
        monkeypatch.undo()

    assert sorted(results) == [2, 4, 6]
    assert completed == [1, 2, 3]
    assert [future._value for future in executor.submitted] == [1, 2, 3]


def test_run_process_pool_propagates_worker_errors() -> None:
    executor = FakeProcessPoolExecutor(max_workers=2)

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(
        workers.concurrent.futures,
        "ProcessPoolExecutor",
        lambda max_workers: executor,
    )
    monkeypatch.setattr(
        workers.concurrent.futures,
        "as_completed",
        lambda futures: iter(futures),
    )

    with pytest.raises(ValueError, match="worker boom"):
        try:
            run_process_pool(fail_on_two, [1, 2], max_workers=2)
        finally:
            monkeypatch.undo()
