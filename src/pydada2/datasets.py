"""Bundled example data (the ``inst/extdata`` files and data sets of the R package)."""

from __future__ import annotations

from functools import cache
from importlib.resources import files
from pathlib import Path

import numpy as np
import numpy.typing as npt

__all__ = [
    "err_balanced_f",
    "err_balanced_r",
    "example_path",
    "list_examples",
    "phix_genome",
    "tperr1",
]


def _data_dir() -> Path:
    return Path(str(files("pydada2") / "data"))


def list_examples() -> list[str]:
    """List the names of the bundled example files.

    Returns
    -------
    list of str
        File names accepted by :func:`example_path`.

    Examples
    --------
    >>> "sam1F.fastq.gz" in list_examples()
    True
    """
    return sorted(p.name for p in (_data_dir() / "extdata").iterdir() if p.is_file())


def example_path(name: str) -> Path:
    """Return the path of a bundled example file.

    Parameters
    ----------
    name : str
        File name, e.g. ``"sam1F.fastq.gz"`` (see :func:`list_examples`).

    Returns
    -------
    pathlib.Path
        Absolute path of the file.

    Raises
    ------
    FileNotFoundError
        If no such example file exists.

    Examples
    --------
    >>> example_path("sam1F.fastq.gz").name
    'sam1F.fastq.gz'
    """
    path = _data_dir() / "extdata" / name
    if not path.is_file():
        raise FileNotFoundError(
            f"No example file named {name!r}. Available: {', '.join(list_examples())}."
        )
    return path


def _load_err(name: str) -> npt.NDArray[np.float64]:
    arr: npt.NDArray[np.float64] = np.load(_data_dir() / f"{name}.npy")
    arr.setflags(write=False)
    return arr


def err_balanced_f() -> npt.NDArray[np.float64]:
    """Return the forward-read error matrix shipped with DADA2 (R: ``errBalancedF``).

    Returns
    -------
    numpy.ndarray
        Read-only array of shape (16, 41): transitions (A2A, A2C, ..., T2T) by quality 0..40.

    Examples
    --------
    >>> err_balanced_f().shape
    (16, 41)
    """
    return _load_err("errBalancedF")


def err_balanced_r() -> npt.NDArray[np.float64]:
    """Return the reverse-read error matrix shipped with DADA2 (R: ``errBalancedR``).

    Returns
    -------
    numpy.ndarray
        Read-only array of shape (16, 41).
    """
    return _load_err("errBalancedR")


def tperr1() -> npt.NDArray[np.float64]:
    """Return the example error matrix ``tperr1`` from the R package.

    Returns
    -------
    numpy.ndarray
        Read-only array of shape (16, 41).
    """
    return _load_err("tperr1")


@cache
def phix_genome() -> str:
    """Return the phiX174 genome sequence used to detect phiX contamination.

    Returns
    -------
    str
        Upper-case genome sequence.

    Examples
    --------
    >>> len(phix_genome())
    5386
    """
    lines = example_path("phix_genome.fa").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()
