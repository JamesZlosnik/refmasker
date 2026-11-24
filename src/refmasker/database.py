"""The geNomad database: download, checksum manifest and verification.

The default prophage detector needs geNomad's database, which is downloaded once and then
reused by every run. refmasker does not fetch it itself: `genomad download-database` does,
so the database always matches the installed geNomad version, and refmasker never has to
hold a download URL that goes stale.

What refmasker adds is a record of exactly what arrived. Every file in the database folder
gets a SHA-256, written next to the folder as a `sha256sum`-compatible manifest, plus one
short *fingerprint* over the whole set. That makes a database checkable later (`--verify`),
and comparable between machines or between everyone sharing one read-only copy on a cluster.

It is not a check against a checksum published by geNomad: none is published beside the
archive, and geNomad's downloader verifies none, so HTTPS to its host is what protects the
transfer. See the README.
"""

from __future__ import annotations

import hashlib
import logging
import shutil
import sys
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from refmasker import __version__, external
from refmasker.safeio import open_replacing

log = logging.getLogger(__name__)

# The folder `genomad download-database DEST` creates inside DEST.
DB_DIRNAME = "genomad_db"

# Where geNomad hosts its database, shown to the user so the download is not a black box.
# refmasker never fetches from here: geNomad's downloader holds the real URL (it reads
# releases.txt for the database version matching the installed geNomad, then fetches
# genomad_db_v<version>.tar.gz) and prints the archive URL it uses, which is the one to
# trust if this ever goes stale. Checked against geNomad 1.13.0.
DB_SOURCE_URL = "https://portal.nersc.gov/genomad/__data__/"

# The manifest is written next to the database folder, not inside it, so the database itself
# can be made read-only (chmod -R a-w) without stopping a later re-download or re-verify.
MANIFEST_SUFFIX = ".sha256sums"

# Download size quoted to the user. From geNomad's documentation; `download` reports the real
# size afterwards from the files on disk.
APPROX_DOWNLOAD_GB = 1.4

# Free space below this only warns: it is a margin for the compressed archive plus its
# extracted contents, and how much the extracted database needs depends on its version.
RECOMMENDED_FREE_GB = 6.0

# Read files in 1 MiB blocks so a multi-GB database is never held in memory.
HASH_BLOCK = 1 << 20


@dataclass(frozen=True)
class FileSum:
    """The SHA-256 and size of one file, named relative to the database folder."""

    path: str  # relative, POSIX separators, so manifests match across platforms
    size: int
    sha256: str


@dataclass
class Manifest:
    """A parsed checksum manifest: what was recorded when the database was downloaded."""

    sums: dict[str, str] = field(default_factory=dict)  # relative path -> sha256
    fingerprint: str = ""
    db_version: str = ""
    written_by: str = ""
    written: str = ""


@dataclass
class VerifyResult:
    """The outcome of checking a database against its manifest."""

    checked: int
    total_bytes: int
    fingerprint: str
    recorded_fingerprint: str
    changed: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    added: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True only if every recorded file is present, unchanged, and nothing was added."""
        return not (self.changed or self.missing or self.added)


def human_bytes(n: int) -> str:
    """Format a byte count for people: '1.4 GB', '612.0 MB'."""
    for unit, size in (("GB", 1 << 30), ("MB", 1 << 20), ("kB", 1 << 10)):
        if n >= size:
            return f"{n / size:.1f} {unit}"
    return f"{n} B"


def manifest_path(db_dir: Path) -> Path:
    """Where the manifest for `db_dir` lives: a sibling file, e.g. genomad_db.sha256sums."""
    return db_dir.with_name(db_dir.name + MANIFEST_SUFFIX)


def resolve_db_dir(path: Path) -> Path:
    """Find the database folder from either the parent folder or the folder itself.

    Accepts what the user is likely to have to hand: the DEST they downloaded into (which
    holds `genomad_db`), or `genomad_db` itself (what REFMASKER_GENOMAD_DB points at).
    """
    inner = path / DB_DIRNAME
    if inner.is_dir():
        return inner
    if path.is_dir() and path.name == DB_DIRNAME:
        return path
    raise FileNotFoundError(
        f"No geNomad database found at {path}: expected {inner} (from `refmasker download-db "
        f"{path}`), or a folder named {DB_DIRNAME}."
    )


def free_space_gb(path: Path) -> float:
    """Free space in GB on the filesystem holding `path`, or the nearest existing parent."""
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe).free / (1 << 30)


class Progress:
    """A progress bar on stderr for a long local read, with a log line when not a terminal.

    Checksumming a 1.4 GB database takes long enough to look like a hang, so it reports as it
    goes: a redrawn bar on an interactive terminal, and progress lines in the log otherwise
    (inside a batch job, where a redrawn bar would fill the job output with partial lines).
    """

    WIDTH = 24

    def __init__(self, total: int, label: str = "", enabled: bool | None = None):
        self.total = total
        self.label = label
        self.done_bytes = 0
        self.enabled = sys.stderr.isatty() if enabled is None else enabled
        self._logged_pct = 0

    def advance(self, n: int) -> None:
        """Count `n` more bytes done and redraw (or log) if there is anything new to show."""
        self.done_bytes += n
        pct = int(100 * self.done_bytes / self.total) if self.total else 100
        if self.enabled:
            filled = pct * self.WIDTH // 100
            bar = "#" * filled + "-" * (self.WIDTH - filled)
            done, total = human_bytes(self.done_bytes), human_bytes(self.total)
            sys.stderr.write(f"\r  [{bar}] {pct:3d}%  {done} / {total}")
            sys.stderr.flush()
        elif pct >= self._logged_pct + 25:
            self._logged_pct = pct - pct % 25
            sys.stdout.flush()  # keep these lines in order with the caller's own output
            log.info("%s %d%% (%s)", self.label or "progress", pct, human_bytes(self.done_bytes))

    def done(self) -> None:
        """Finish the line, so later output does not land on top of the bar."""
        if self.enabled:
            sys.stderr.write("\n")
            sys.stderr.flush()


def db_files(db_dir: Path) -> list[Path]:
    """Every regular file in the database folder, in a stable order."""
    return sorted(p for p in db_dir.rglob("*") if p.is_file())


def hash_file(path: Path, bar: Progress | None = None) -> tuple[str, int]:
    """SHA-256 one file, in blocks so a multi-GB file is never held in memory.

    Returns the digest and the number of bytes read, and counts those bytes towards `bar`.
    """
    h = hashlib.sha256()
    size = 0
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(HASH_BLOCK), b""):
            h.update(block)
            size += len(block)
            if bar is not None:
                bar.advance(len(block))
    return h.hexdigest(), size


def hash_tree(db_dir: Path, show_progress: bool | None = None) -> list[FileSum]:
    """SHA-256 every file in the database folder, reporting progress as it reads."""
    files = db_files(db_dir)
    if not files:
        raise FileNotFoundError(f"{db_dir} contains no files")
    total = sum(p.stat().st_size for p in files)
    bar = Progress(total, label="checksumming", enabled=show_progress)
    sums = []
    for path in files:
        digest, size = hash_file(path, bar)
        sums.append(FileSum(path.relative_to(db_dir).as_posix(), size, digest))
    bar.done()
    return sums


def fingerprint(sums: Iterable[FileSum]) -> str:
    """One SHA-256 over every file's checksum and name: a short name for a whole database.

    Two databases with the same fingerprint hold byte-identical files under the same names,
    so two machines (or two users of a shared copy) can compare one 64-character string
    instead of thousands of checksums.
    """
    h = hashlib.sha256()
    for fs in sorted(sums, key=lambda f: f.path):
        h.update(f"{fs.sha256}  {fs.path}\n".encode())
    return h.hexdigest()


def write_manifest(db_dir: Path, sums: list[FileSum], db_version: str) -> Path:
    """Write the checksum manifest beside the database folder and return its path.

    The checksum lines are the `sha256sum` format ('<hex>  <path>', two spaces, paths
    relative to the database folder), so GNU coreutils can check them too. The header lines
    are comments, which `sha256sum -c` skips.
    """
    path = manifest_path(db_dir)
    total = sum(f.size for f in sums)
    fp = fingerprint(sums)
    written = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    with open_replacing(path) as fh:
        fh.write(f"# refmasker {__version__} checksum manifest for a geNomad database\n")
        fh.write(f"# database folder: {db_dir.name}\n")
        fh.write(f"# geNomad database version: {db_version}\n")
        fh.write(f"# files: {len(sums)}  total bytes: {total}\n")
        fh.write(f"# written: {written}\n")
        fh.write(f"# fingerprint (SHA-256 of the lines below): {fp}\n")
        fh.write("#\n")
        fh.write(f"# Re-check with: refmasker download-db {db_dir.parent} --verify\n")
        fh.write(f"# or: cd {db_dir} && sha256sum -c ../{path.name}\n")
        for fs in sorted(sums, key=lambda f: f.path):
            fh.write(f"{fs.sha256}  {fs.path}\n")
    return path


def read_manifest(path: Path) -> Manifest:
    """Parse a manifest written by write_manifest (or any sha256sum-format file)."""
    m = Manifest()
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            # Header comments: pick out the values refmasker records.
            text = line.lstrip("# ")
            low = text.lower()
            if low.startswith("fingerprint") and ":" in text:
                m.fingerprint = text.split(":", 1)[1].strip()
            elif low.startswith("genomad database version:"):
                m.db_version = text.split(":", 1)[1].strip()
            elif low.startswith("refmasker ") and "manifest" in low:
                m.written_by = text.split()[1]
            elif low.startswith("written:"):
                m.written = text.split(":", 1)[1].strip()
            continue
        if "  " not in line:
            raise ValueError(f"{path}:{lineno}: not a '<sha256>  <path>' line")
        digest, name = line.split("  ", 1)
        m.sums[name.strip()] = digest.strip()
    if not m.sums:
        raise ValueError(f"{path}: no checksum lines found")
    return m


def verify(db_dir: Path, show_progress: bool | None = None) -> tuple[VerifyResult, Manifest]:
    """Re-checksum a database and compare it with its manifest.

    Reports files whose contents changed, files that are gone, and files that appeared, and
    recomputes the fingerprint so it can be compared with the recorded one.
    """
    path = manifest_path(db_dir)
    if not path.is_file():
        raise FileNotFoundError(
            f"No checksum manifest at {path}, so this database cannot be verified. It was "
            "either downloaded by a refmasker version that did not write one, or by geNomad "
            "directly. Re-run `refmasker download-db` into a new folder to get a manifest."
        )
    recorded = read_manifest(path)
    sums = hash_tree(db_dir, show_progress=show_progress)
    found = {f.path: f.sha256 for f in sums}
    changed = sorted(p for p, h in found.items() if p in recorded.sums and h != recorded.sums[p])
    missing = sorted(set(recorded.sums) - set(found))
    added = sorted(set(found) - set(recorded.sums))
    result = VerifyResult(
        checked=len(sums),
        total_bytes=sum(f.size for f in sums),
        fingerprint=fingerprint(sums),
        recorded_fingerprint=recorded.fingerprint,
        changed=changed,
        missing=missing,
        added=added,
    )
    return result, recorded


def download(dest: Path) -> Path:
    """Run `genomad download-database DEST` and return the database folder it created.

    geNomad's own output goes straight to the terminal (see external.run_attached), so the
    user sees which URL it fetches and whatever progress it reports, rather than waiting on
    a silent command for a 1.4 GB download.
    """
    external.run_attached(["genomad", "download-database", str(dest)])
    db_dir = dest / DB_DIRNAME
    if not db_dir.is_dir():
        raise external.ToolError(
            f"geNomad's downloader finished but there is no {db_dir}. Check its output above."
        )
    return db_dir
