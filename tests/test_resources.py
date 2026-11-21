import logging
import os

import pytest

from refmasker.resources import limit_library_threads, resolve_threads


def test_default_is_one_thread():
    assert resolve_threads(None, {}) == 1


def test_scheduler_allocation():
    assert resolve_threads(None, {"SLURM_CPUS_PER_TASK": "8"}) == 8
    assert resolve_threads(None, {"NSLOTS": "4"}) == 4
    assert resolve_threads(None, {"SLURM_CPUS_ON_NODE": "6"}) == 6  # no --cpus-per-task
    assert resolve_threads(None, {"NCPUS": "3"}) == 3  # PBS Pro
    assert resolve_threads(None, {"SLURM_CPUS_PER_TASK": "2", "SLURM_CPUS_ON_NODE": "8"}) == 2


def test_explicit_wins():
    assert resolve_threads(2, {"SLURM_CPUS_PER_TASK": "8"}) == 2


def test_warns_when_asking_for_more_than_the_job_has(caplog):
    with caplog.at_level(logging.WARNING):
        assert resolve_threads(16, {"SLURM_CPUS_PER_TASK": "4"}) == 16
    assert "this job was given" in caplog.text


def test_limit_library_threads(monkeypatch):
    env = {"JAVA_TOOL_OPTIONS": "-Xmx2g"}
    monkeypatch.setattr(os, "environ", env)
    limit_library_threads(3)
    assert env["OMP_NUM_THREADS"] == "3"
    assert env["TF_NUM_INTRAOP_THREADS"] == "3"
    assert env["TF_NUM_INTEROP_THREADS"] == "1"
    assert env["JAVA_TOOL_OPTIONS"] == "-Xmx2g -XX:ActiveProcessorCount=3"
    limit_library_threads(3)  # running twice does not add the option twice
    assert env["JAVA_TOOL_OPTIONS"].count("ActiveProcessorCount") == 1


def test_rejects_zero():
    with pytest.raises(ValueError):
        resolve_threads(0, {})
