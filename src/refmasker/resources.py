"""Decide how many threads to use without hogging a shared machine."""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping

log = logging.getLogger(__name__)

# Environment variables set by schedulers to the number of cores a job was given, in the
# order they are checked: SLURM (--cpus-per-task, then the CPUs allocated on this node, which
# is set even without --cpus-per-task), SGE, PBS Pro, Torque.
SCHEDULER_VARS = ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "NSLOTS", "NCPUS", "PBS_NUM_PPN")

# Variables that stop numerical libraries (OpenMP, BLAS) spawning one thread per core.
LIBRARY_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def scheduler_allocation(environ: Mapping[str, str]) -> tuple[str, int] | None:
    """The (variable, cores) of the first scheduler variable that is set, or None outside a job."""
    for var in SCHEDULER_VARS:
        val = environ.get(var, "").strip()
        if val.isdigit() and int(val) > 0:
            return var, int(val)
    return None


def resolve_threads(requested: int | None, environ: Mapping[str, str] | None = None) -> int:
    """Return the thread count to use.

    An explicit --threads wins. Otherwise use the scheduler allocation if there is one,
    else 1. Never default to every core on the machine.
    """
    env = os.environ if environ is None else environ
    allocation = scheduler_allocation(env)
    # An explicit request is used as given, after a sanity check and a warning if it is more
    # than the job was given (or, outside a job, more than the machine has).
    if requested is not None:
        if requested < 1:
            raise ValueError("--threads must be at least 1")
        ncpu = os.cpu_count() or 1
        if allocation is not None and requested > allocation[1]:
            log.warning(
                "--threads %d is more than the %d CPUs this job was given (%s); the job will "
                "compete with itself and, without CPU limits on the node, with other users",
                requested,
                allocation[1],
                allocation[0],
            )
        elif requested > ncpu:
            log.warning("--threads %d is more than the %d CPUs on this machine", requested, ncpu)
        return requested
    # No request: use the job's allocation.
    if allocation is not None:
        log.info("Using %s=%d threads from the job scheduler", *allocation)
        return allocation[1]
    # Not inside a scheduler job: use one thread to be safe on a shared machine.
    return 1


def limit_library_threads(threads: int) -> None:
    """Cap thread pools of numerical libraries in this process and its children.

    Child processes (the external tools) inherit os.environ, so the caps apply to them too.
    """
    for var in LIBRARY_THREAD_VARS:
        os.environ[var] = str(threads)
    # TensorFlow (used by geNomad's neural-network step): one operation at a time, each using
    # at most `threads` threads, so the total stays near `threads`.
    os.environ["TF_NUM_INTRAOP_THREADS"] = str(threads)
    os.environ["TF_NUM_INTEROP_THREADS"] = "1"
    # Java (MinCED) sizes its garbage-collector and other thread pools from the number of CPUs
    # it sees; make it see `threads`. Keep any options already set, and a user's own setting.
    java_opts = os.environ.get("JAVA_TOOL_OPTIONS", "")
    if "ActiveProcessorCount" not in java_opts:
        os.environ["JAVA_TOOL_OPTIONS"] = f"{java_opts} -XX:ActiveProcessorCount={threads}".strip()
