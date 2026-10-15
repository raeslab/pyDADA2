"""``n_jobs`` handling."""

from __future__ import annotations

import os

__all__ = ["effective_n_jobs"]


def effective_n_jobs(n_jobs: int | None) -> int:
    """Translate an ``n_jobs`` argument into a positive number of threads.

    Parameters
    ----------
    n_jobs : int or None
        ``None`` means 1. Negative values count back from the number of CPUs
        (``-1`` = all CPUs, ``-2`` = all but one), as in scikit-learn.

    Returns
    -------
    int
        Number of threads to use (at least 1).

    Raises
    ------
    ValueError
        If ``n_jobs`` is 0 or not an integer.
    """
    if n_jobs is None:
        return 1
    if isinstance(n_jobs, bool) or not isinstance(n_jobs, int):
        raise ValueError(f"n_jobs must be an int or None, got {n_jobs!r}.")
    if n_jobs == 0:
        raise ValueError("n_jobs == 0 has no meaning; use None, a positive int or -1.")
    if n_jobs > 0:
        return n_jobs
    n_cpu = _cpu_count()
    return max(1, n_cpu + 1 + n_jobs)


def _cpu_count() -> int:
    n: int | None
    if hasattr(os, "process_cpu_count"):  # Python >= 3.13
        n = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):
        n = len(os.sched_getaffinity(0))
    else:  # pragma: no cover - macOS/Windows on Python < 3.13
        n = os.cpu_count()
    return n or 1
