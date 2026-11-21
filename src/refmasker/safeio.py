"""Writing output files safely in folders that other users might be able to reach."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any


def current_umask() -> int:
    """The process umask. os.umask can only be read by setting it, so set it straight back."""
    mask = os.umask(0)
    os.umask(mask)
    return mask


@contextmanager
def open_replacing(
    path: str | Path, mode: str = "w", newline: str | None = None
) -> Iterator[IO[Any]]:
    """Open a temporary file next to `path` for writing; move it into place when done.

    os.replace swaps the directory entry itself, so if `path` is a symbolic link (for example
    one planted by another user) the link is replaced rather than followed, and nobody ever
    sees a half-written file. If writing fails, the temporary file is removed and `path` is
    left as it was. The file gets the same permissions as with open(): 0o666 minus the umask.
    """
    path = Path(path)
    # mkstemp creates the file exclusively (it never opens an existing file or link) and
    # returns an open descriptor, so the file is not looked up by name again before writing.
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, mode, newline=newline) as fh:
            os.fchmod(fh.fileno(), 0o666 & ~current_umask())
            yield fh
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
