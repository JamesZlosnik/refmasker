import pytest

from refmasker.resources import resolve_threads


def test_default_is_one_thread():
    assert resolve_threads(None, {}) == 1


def test_scheduler_allocation():
    assert resolve_threads(None, {"SLURM_CPUS_PER_TASK": "8"}) == 8
    assert resolve_threads(None, {"NSLOTS": "4"}) == 4


def test_explicit_wins():
    assert resolve_threads(2, {"SLURM_CPUS_PER_TASK": "8"}) == 2


def test_rejects_zero():
    with pytest.raises(ValueError):
        resolve_threads(0, {})
