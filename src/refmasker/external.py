"""Running and checking external command-line tools."""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

# How to ask each tool for its version.
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
    pass


def missing_tools(tools: list[str]) -> list[str]:
    return [t for t in tools if shutil.which(t) is None]


def tool_version(tool: str) -> str:
    cmd = VERSION_COMMANDS.get(tool)
    if not cmd or shutil.which(cmd[0]) is None:
        return "installed" if shutil.which(tool) else "not found"
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
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


def run(
    cmd: list[str],
    log_dir: Path,
    name: str,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess:
    """Run a command, logging it and saving stderr to `log_dir/name.log`.

    If `stdout_path` is given, stdout is written there; otherwise it is captured.
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
        tail = "".join(err_path.read_text().splitlines(keepends=True)[-20:])
        raise ToolError(
            f"{cmd[0]} failed with exit code {p.returncode}. Last lines of {err_path}:\n{tail}"
        )
    return p
