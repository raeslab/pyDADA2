import pytest

from pydada2._parallel import _cpu_count, effective_n_jobs


def test_effective_n_jobs():
    n = _cpu_count()
    assert effective_n_jobs(None) == 1
    assert effective_n_jobs(3) == 3
    assert effective_n_jobs(-1) == n
    assert effective_n_jobs(-2) == max(1, n - 1)
    assert effective_n_jobs(-10_000) == 1


@pytest.mark.parametrize("bad", [0, 1.5, "2", True])
def test_effective_n_jobs_invalid(bad):
    with pytest.raises(ValueError):
        effective_n_jobs(bad)
