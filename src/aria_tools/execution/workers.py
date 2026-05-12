"""Native Python worker execution helpers."""

from __future__ import annotations

import concurrent.futures
import os
from collections.abc import Callable, Iterable
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def resolve_worker_count(num_workers: int | str) -> int:
    """Return a concrete worker count from CLI-style worker settings."""

    if isinstance(num_workers, str) and num_workers.upper() == "ALL_CPUS":
        return os.cpu_count() or 1
    return int(num_workers)


def run_process_pool(
    worker: Callable[[T], R],
    work_items: Iterable[T],
    *,
    max_workers: int,
    on_complete: Callable[[int], None] | None = None,
) -> list[R]:
    """Run work items in a process pool and collect results in memory."""

    results: list[R] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(worker, item) for item in work_items]
        for completed, future in enumerate(
            concurrent.futures.as_completed(futures), start=1
        ):
            results.append(future.result())
            if on_complete is not None:
                on_complete(completed)

    return results
