"""Decide how many threads to use without hogging a shared machine."""

from __future__ import annotations

import logging
import os

log = logging.getLogger(__name__)

# Environment variables set by schedulers to the number of cores a job was given.
SCHEDULER_VARS = ("SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NUM_PPN")

# Variables that stop numerical libraries (OpenMP, BLAS) spawning one thread per core.
LIBRARY_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def resolve_threads(requested: int | None, environ: dict[str, str] | None = None) -> int:
    """Return the thread count to use.

    An explicit --threads wins. Otherwise use the scheduler allocation if there is one,
    else 1. Never default to every core on the machine.
    """
    env = os.environ if environ is None else environ
    # An explicit request is used as given (after a sanity check and a warning if too high).
    if requested is not None:
        if requested < 1:
            raise ValueError("--threads must be at least 1")
        ncpu = os.cpu_count() or 1
        if requested > ncpu:
            log.warning("--threads %d is more than the %d CPUs on this machine", requested, ncpu)
        return requested
    # No request: use the job's allocation, checking SLURM, then SGE, then PBS.
    for var in SCHEDULER_VARS:
        val = env.get(var, "").strip()
        if val.isdigit() and int(val) > 0:
            log.info("Using %s=%s threads from the job scheduler", var, val)
            return int(val)
    # Not inside a scheduler job: use one thread to be safe on a shared machine.
    return 1


def limit_library_threads(threads: int) -> None:
    """Cap thread pools of numerical libraries in this process and its children."""
    for var in LIBRARY_THREAD_VARS:
        os.environ[var] = str(threads)
