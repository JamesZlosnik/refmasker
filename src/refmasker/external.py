"""Running and checking external command-line tools."""

from __future__ import annotations

import logging
import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

log = logging.getLogger(__name__)

# How to ask each tool for its version. None means the tool has no version option and is just
# reported as "installed" (show-coords comes with MUMmer, so nucmer's version covers it).
VERSION_COMMANDS = {
    "nucmer": ["nucmer", "--version"],
    "show-coords": None,
    "isescan.py": ["isescan.py", "--version"],
    "genomad": ["genomad", "--version"],
    "PhiSpy.py": ["PhiSpy.py", "--version"],
    "barrnap": ["barrnap", "--version"],
    "dustmasker": ["dustmasker", "-version-full"],
    "minced": ["minced", "--version"],
}


class ToolError(RuntimeError):
    """An external tool is missing, failed, or produced no usable output."""


def missing_tools(tools: list[str]) -> list[str]:
    """Return the tools in `tools` that are not on PATH."""
    return [t for t in tools if shutil.which(t) is None]


def writable_by_others(path: str | Path) -> bool:
    """True if any user could change `path` or the folder holding it (e.g. swap a program).

    World-writable folders with the sticky bit (like /tmp) are not counted: there, users
    cannot rename or delete each other's files.
    """
    p = Path(path)
    for item in (p, p.parent):
        try:
            mode = os.stat(item).st_mode
        except OSError:
            continue
        sticky_dir = stat.S_ISDIR(mode) and mode & stat.S_ISVTX
        if mode & stat.S_IWOTH and not sticky_dir:
            return True
    return False


def tool_version(tool: str) -> str:
    """Return a tool's version number, or "not found", "installed" or "unknown".

    The version is the first number like X.Y or X.Y.Z in the tool's output.
    """
    cmd = VERSION_COMMANDS.get(tool)
    # No known version command, or its program is missing: just say whether the tool exists.
    if not cmd or shutil.which(cmd[0]) is None:
        return "installed" if shutil.which(tool) else "not found"
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    # Some tools print their version on stderr, so search both streams.
    text = (p.stdout or "") + "\n" + (p.stderr or "")
    for line in text.splitlines():
        line = line.strip()
        # Skip JVM chatter (minced is Java) and blank lines.
        if not line or line.startswith(("Picked up", "WARNING")):
            continue
        m = re.search(r"\d+\.\d+(\.\d+)?", line)
        if m:
            return m.group(0)
    return "unknown"


def run_attached(cmd: list[str], cwd: Path | None = None) -> None:
    """Run a command with its output going straight to this terminal, and wait for it.

    Unlike `run`, nothing is captured or redirected: the child inherits stdout and stderr, so
    a tool that draws a progress bar or prints the URL it is fetching shows that to the user
    as it happens. Used for the long, interactive geNomad database download, where a silent
    command looks like a hang. Nothing is logged to a file, so the output lives only in the
    terminal (or the job's output file).
    """
    # Flush first: our own stdout may be block-buffered (a job log), and anything still in
    # the buffer would otherwise appear after this log line and after the child's output.
    sys.stdout.flush()
    sys.stderr.flush()
    log.info("Running: %s", " ".join(str(c) for c in cmd))
    sys.stderr.flush()
    try:
        p = subprocess.run(cmd, cwd=cwd)
    except OSError as exc:
        raise ToolError(f"Could not run {cmd[0]}: {exc}") from exc
    if p.returncode != 0:
        raise ToolError(f"{cmd[0]} failed with exit code {p.returncode}; its output is above.")


def run(
    cmd: list[str],
    log_dir: Path,
    name: str,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess:
    """Run a command, logging it and saving stderr to `log_dir/name.log`.

    If `stdout_path` is given, stdout is written there; otherwise it is captured.
    The command is an argument list run without a shell, so spaces or special characters in
    file names cannot be misread as shell syntax.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    err_path = log_dir / f"{name}.log"
    log.info("Running: %s", " ".join(str(c) for c in cmd))
    with open(err_path, "w") as err:
        if stdout_path is not None:
            with open(stdout_path, "w") as out:
                p = subprocess.run(cmd, stdout=out, stderr=err, text=True, cwd=cwd)
        else:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=err, text=True, cwd=cwd)
    if p.returncode != 0:
        # Put the end of the tool's log in the error message so the user can see why it failed.
        tail = "".join(err_path.read_text().splitlines(keepends=True)[-20:])
        raise ToolError(
            f"{cmd[0]} failed with exit code {p.returncode}. Last lines of {err_path}:\n{tail}"
        )
    return p
