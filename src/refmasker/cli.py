"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import time
from pathlib import Path

from refmasker import __version__

log = logging.getLogger(__name__)

# A read longer than a repeat can be placed by the unique sequence on either side, so by default
# only repeats at least as long as the reads are masked, and never shorter than 200 bp.
DEFAULT_READ_LENGTH = 150
MIN_REPEAT_FLOOR = 200

# Example commands shown at the end of `refmasker --help`.
EPILOG = """\
examples:
  refmasker run ref.fna -o ref_mask
  refmasker run ref.fna --annotation ref.gff -o ref_mask --threads 4
  refmasker run ref.fna -o ref_mask --categories default,crispr,lowcomplexity
  refmasker run ref.fna -o ref_mask --read-length 300 --formats html
  refmasker compare ref_mask/mask.bed published_mask.bed
  refmasker check
  refmasker download-db /path/to/dbs            (the geNomad database, once)
  refmasker download-db /path/to/dbs --verify   (re-check its checksums)
"""


def _positive_int(text: str) -> int:
    """argparse type: an integer >= 0 (zero is allowed, despite the name)."""
    v = int(text)
    if v < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return v


def _identity(text: str) -> float:
    """argparse type: a percent identity, above 0 and at most 100."""
    v = float(text)
    if not 0 < v <= 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return v


def build_parser() -> argparse.ArgumentParser:
    """Define the command line: the `run`, `compare`, `download-db` and `check` subcommands."""
    # The category names and defaults in the help text come from the detector registry.
    from refmasker.detectors import categories, default_categories

    p = argparse.ArgumentParser(
        prog="refmasker",
        description="Build a BED mask of problematic regions in a bacterial reference genome "
        "for snippy core-genome SNP analysis.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"refmasker {__version__}")
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    # "repeats (default); is (default); ...", listed in the --categories help.
    cat_help = "; ".join(
        f"{c.name}{' (default)' if c.default else ''}" for c in categories().values()
    )

    # ---- refmasker run ----
    r = sub.add_parser(
        "run",
        help="build a mask for one reference genome",
        description="Build a mask for one reference genome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    r.add_argument("reference", help="reference genome FASTA (may be gzipped)")
    r.add_argument("-o", "--outdir", required=True, help="output folder")
    r.add_argument("--annotation", help="GFF3 or GenBank annotation for the same reference")
    r.add_argument(
        "--categories",
        default=",".join(default_categories()),
        help=f"comma-separated categories to mask; also 'default' and 'all'. Available: {cat_help}",
    )
    g = r.add_argument_group("repeats")
    g.add_argument(
        "--read-length",
        type=int,
        default=DEFAULT_READ_LENGTH,
        help="read length the mask is built for; sets the default --repeat-min-len",
    )
    g.add_argument(
        "--repeat-min-len",
        type=_positive_int,
        default=None,
        help=f"minimum repeat length in bp [default: max({MIN_REPEAT_FLOOR}, --read-length)]",
    )
    # "%%" is how argparse help text shows a literal "%".
    g.add_argument(
        "--repeat-min-id", type=_identity, default=90.0, help="minimum %% identity between copies"
    )
    g = r.add_argument_group("other categories")
    g.add_argument(
        "--phage-tool",
        choices=["genomad", "phispy"],
        default="genomad",
        help="prophage caller (PhiSpy needs a GenBank --annotation)",
    )
    g.add_argument("--genomad-db", help="geNomad database folder [default: $REFMASKER_GENOMAD_DB]")
    g.add_argument(
        "--isescan-partial",
        action="store_true",
        help="also mask IS elements ISEScan calls partial (more remnants, more false positives)",
    )
    g.add_argument(
        "--mask-genes-matching",
        metavar="REGEX",
        help="mask annotated genes whose product, gene name or locus tag matches this regular "
        "expression, e.g. 'PE_PGRS|PPE' (adds category 'genes'; needs --annotation)",
    )
    g.add_argument(
        "--rrna-operon-gap",
        type=_positive_int,
        default=1000,
        help="join rRNA genes this close into one operon region (0 keeps genes separate)",
    )
    g.add_argument(
        "--homopolymer-min-len",
        type=int,
        default=8,
        help="shortest homopolymer run masked by the lowcomplexity category",
    )
    g.add_argument(
        "--extra-bed",
        action="append",
        metavar="BED",
        help="extra regions to mask (category 'user'); can be repeated",
    )
    g.add_argument(
        "--mask-plasmids",
        action="store_true",
        help="mask whole plasmid records (adds category 'plasmid')",
    )
    g.add_argument(
        "--plasmid-ids",
        type=lambda s: [x for x in s.split(",") if x],
        default=[],
        help="comma-separated record IDs to treat as plasmids, in addition to records whose "
        "FASTA header mentions 'plasmid'",
    )
    g = r.add_argument_group("output and resources")
    g.add_argument(
        "--pad", type=_positive_int, default=0, help="bp added to each side of every region"
    )
    g.add_argument(
        "--formats",
        default="html,pdf,png",
        help="report formats: any of html, pdf, png, or 'none'",
    )
    g.add_argument(
        "--threads",
        type=int,
        default=None,
        help="threads to use [default: the job's allocation under SLURM/SGE/PBS, otherwise 1]",
    )
    g.add_argument(
        "--max-genome-mb",
        type=_positive_int,
        default=50,
        help="refuse a reference with more sequence than this (Mb), so a wrong input such as a "
        "metagenome assembly fails at once instead of using a lot of memory and CPU; 0 = no limit",
    )
    g.add_argument(
        "--keep-work", action="store_true", help="keep intermediate tool outputs in OUTDIR/work"
    )
    g.add_argument(
        "--force",
        action="store_true",
        help="overwrite a previous refmasker output folder (refuses folders refmasker did not create)",
    )
    g.add_argument(
        "--redact",
        action="store_true",
        help="for sharing: reports (HTML, PDF, run_info.json) show input file names instead of "
        "full paths and omit FASTA header descriptions. Record IDs are kept (mask.bed needs "
        "them) and refmasker.log keeps full paths",
    )
    g.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")

    # ---- refmasker compare ----
    c = sub.add_parser(
        "compare",
        help="compare two BED masks",
        description="Compare two BED masks: bp in common, Jaccard index, regions unique to each.",
    )
    c.add_argument("a", help="first BED (e.g. refmasker's mask.bed)")
    c.add_argument("b", help="second BED (e.g. a published mask)")
    c.add_argument(
        "--rename",
        action="append",
        default=[],
        metavar="OLD=NEW",
        help="rename a sequence in both files before comparing; can be repeated",
    )
    c.add_argument("--out-prefix", help="write PREFIX.only_a.bed and PREFIX.only_b.bed")

    # ---- refmasker download-db ----
    d = sub.add_parser(
        "download-db",
        help="download the geNomad database (~1.4 GB), or verify a downloaded one",
        description="Download the geNomad database used by the prophage category, then "
        "checksum every file so the database can be verified later.",
    )
    d.add_argument("dest", help="folder to download into; the database goes in DEST/genomad_db")
    d.add_argument(
        "--verify",
        action="store_true",
        help="do not download: re-checksum the database already in DEST and compare it with "
        "the manifest written when it was downloaded",
    )

    # ---- refmasker check ----
    sub.add_parser("check", help="list the external tools refmasker uses and their versions")
    return p


def setup_logging(quiet: bool = False) -> None:
    """Log to stderr; with --quiet only warnings and errors reach the terminal."""
    logging.basicConfig(
        level=logging.WARNING if quiet else logging.INFO,
        format="[refmasker] %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    # The root logger stays at INFO so the run's log file (added later by the pipeline) still
    # records everything; --quiet only raises the level of the terminal handler.
    logging.getLogger().setLevel(logging.INFO)
    if quiet:
        logging.getLogger().handlers[0].setLevel(logging.WARNING)
    # matplotlib and fontTools log a lot of detail while drawing figures; keep them quiet.
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("fontTools").setLevel(logging.WARNING)


def cmd_run(args) -> int:
    """`refmasker run`: fill in defaults that depend on other options, then run the pipeline."""
    from refmasker.pipeline import run
    from refmasker.report import parse_formats

    # --repeat-min-len defaults to the read length, but never below MIN_REPEAT_FLOOR.
    if args.repeat_min_len is None:
        args.repeat_min_len = max(MIN_REPEAT_FLOOR, args.read_length)
    args.formats = parse_formats(args.formats)
    run(args)
    return 0


def cmd_compare(args) -> int:
    """`refmasker compare`: print overlap statistics for two BED masks."""
    from refmasker.compare import compare_files, format_report
    from refmasker.regions import write_mask_bed

    # Turn the repeated --rename OLD=NEW options into a dict.
    rename = {}
    for item in args.rename:
        if "=" not in item:
            raise ValueError(f"--rename expects OLD=NEW, got '{item}'")
        old, new = item.split("=", 1)
        rename[old] = new
    c = compare_files(Path(args.a), Path(args.b), rename)
    print(format_report(c, Path(args.a).name, Path(args.b).name))
    # Optionally save the regions found in only one of the two masks.
    if args.out_prefix:
        write_mask_bed(f"{args.out_prefix}.only_a.bed", c.only_a)
        write_mask_bed(f"{args.out_prefix}.only_b.bed", c.only_b)
    return 0


# Printed after a download, and worth keeping accurate: refmasker records what it received,
# which is not the same as checking the bytes against a checksum geNomad published.
CHECKSUM_NOTE = """\
The fingerprint is one SHA-256 over every file's checksum and name. It shows this copy is
intact, and lets you check that two machines -- or everyone sharing one install -- hold the
same database: compare fingerprints, or re-run with --verify. It cannot show the download
itself was genuine: geNomad publishes no checksum beside the archive and its downloader
verifies none, so HTTPS to its host is what protects the transfer. Treat this fingerprint as
the baseline for a fresh database and compare against it from now on."""


def _duration(seconds: float) -> str:
    """Format an elapsed time as '42s' or '7m12s'."""
    minutes, secs = divmod(int(seconds), 60)
    return f"{minutes}m{secs:02d}s" if minutes else f"{secs}s"


def cmd_download_db(args) -> int:
    """`refmasker download-db`: download the geNomad database, or verify a downloaded one."""
    from refmasker import database

    dest = Path(args.dest)
    if args.verify:
        return _verify_db(database.resolve_db_dir(dest))
    return _download_db(dest)


def _download_db(dest: Path) -> int:
    """Explain what is about to happen, run geNomad's downloader, then checksum the result."""
    from refmasker import database, external
    from refmasker.detectors.prophages import genomad_db_version

    db_dir = dest / database.DB_DIRNAME
    # Never overwrite or delete an existing database: it is ~1.4 GB and may be shared. Checked
    # before the tool check, so "you already have one" wins over "install geNomad".
    if db_dir.exists():
        raise FileExistsError(
            f"{db_dir} already exists. refmasker will not overwrite a database that others may "
            f"be using: check it with `refmasker download-db {dest} --verify`, download into "
            "another folder, or remove the old one yourself first."
        )
    if external.missing_tools(["genomad"]):
        raise external.ToolError("genomad is not on PATH; install it from environment.yml first")

    free = database.free_space_gb(dest)
    manifest_name = database.manifest_path(db_dir).name
    print(f"refmasker {__version__}: geNomad database download")
    print()
    print(f"  Destination    {db_dir}")
    print("  Needed by      the 'prophage' mask category, whose default caller is geNomad")
    print("  Downloader     `genomad download-database`, geNomad's own, so the database always")
    print("                 matches the installed geNomad version")
    print(f"  geNomad        {external.tool_version('genomad')} at {shutil.which('genomad')}")
    print(f"  Source         {database.DB_SOURCE_URL}")
    print("                 geNomad's own host; it prints the exact archive URL below")
    print(f"  Download size  about {database.APPROX_DOWNLOAD_GB} GB, once (not per run)")
    print(f"  Free space     {free:.1f} GB where the database will go")
    print("  Afterwards     refmasker checksums every file with SHA-256 and writes")
    print(f"                 {manifest_name} next to the database folder")
    print()

    # Fail now if it cannot possibly fit; warn when the margin for extraction looks thin.
    if free < database.APPROX_DOWNLOAD_GB:
        raise ValueError(
            f"Only {free:.1f} GB free where {dest} will be created, and the download alone is "
            f"about {database.APPROX_DOWNLOAD_GB} GB. Choose a folder with more space."
        )
    if free < database.RECOMMENDED_FREE_GB:
        log.warning(
            "%.1f GB free: geNomad downloads a compressed archive and then extracts it, so "
            "about %.0f GB is a safer margin",
            free,
            database.RECOMMENDED_FREE_GB,
        )

    dest.mkdir(parents=True, exist_ok=True)
    print("Downloading. Everything between the lines is geNomad's own output:")
    print("-" * 78)
    t0 = time.time()
    db_dir = database.download(dest)
    print("-" * 78)
    print(f"Download finished in {_duration(time.time() - t0)}.")
    print()

    files = database.db_files(db_dir)
    total = sum(f.stat().st_size for f in files)
    print(
        f"Checksumming {len(files)} files ({database.human_bytes(total)}) with SHA-256. "
        "This reads the whole database once:"
    )
    sums = database.hash_tree(db_dir)
    db_version = genomad_db_version(db_dir)
    manifest = database.write_manifest(db_dir, sums, db_version)
    print()
    print(f"  Database version   {db_version}")
    print(f"  Files              {len(sums)} ({database.human_bytes(total)})")
    print(f"  Fingerprint        {database.fingerprint(sums)}")
    print(f"  Manifest           {manifest}")
    print()
    print(CHECKSUM_NOTE)
    print()
    print("Next:")
    print(f"  export REFMASKER_GENOMAD_DB={db_dir}")
    print("      (add that line to ~/.bashrc, or pass --genomad-db to `refmasker run`)")
    print("  refmasker check")
    print("      confirms refmasker finds the database and shows its version")
    print(f"  refmasker download-db {dest} --verify")
    print("      re-checks every file against the manifest, any time")
    print()
    print("On a shared cluster, make the database read-only so everyone uses the same copy")
    print(f"and nobody can change it:  chmod -R a-w {db_dir}")
    return 0


def _verify_db(db_dir: Path) -> int:
    """Re-checksum a downloaded database and report whether it still matches its manifest."""
    from refmasker import database

    print(f"refmasker {__version__}: verifying a geNomad database")
    print()
    print(f"  Database  {db_dir}")
    print(f"  Manifest  {database.manifest_path(db_dir)}")
    print()
    files = database.db_files(db_dir)
    total = sum(f.stat().st_size for f in files)
    print(
        f"Re-checksumming {len(files)} files ({database.human_bytes(total)}) with SHA-256. "
        "This reads the whole database once:"
    )
    result, recorded = database.verify(db_dir)
    # An old or hand-made manifest may have no fingerprint line; only compare if it has one.
    fingerprint_differs = (
        bool(result.recorded_fingerprint) and result.fingerprint != result.recorded_fingerprint
    )
    print()
    if result.ok and not fingerprint_differs:
        print(f"  OK: {result.checked} files, every checksum matches the manifest.")
        print(f"  Fingerprint       {result.fingerprint}")
        print(f"  Database version  {recorded.db_version or 'unknown'}")
        by = recorded.written_by or "an unknown version"
        print(f"  Recorded by       refmasker {by} on {recorded.written or 'an unknown date'}")
        return 0
    print("  FAILED: this database does not match its manifest.")
    for label, items in (
        ("changed", result.changed),
        ("missing", result.missing),
        ("not in the manifest", result.added),
    ):
        if items:
            shown = ", ".join(items[:10])
            if len(items) > 10:
                shown += f", +{len(items) - 10} more"
            print(f"    {len(items)} file(s) {label}: {shown}")
    if fingerprint_differs:
        print(f"    fingerprint now       {result.fingerprint}")
        print(f"    fingerprint recorded  {result.recorded_fingerprint}")
    print()
    print("  A database that does not match may be incomplete, corrupted, a different version,")
    print("  or changed since it was downloaded. Download a fresh copy into a new folder, and")
    print("  do not use this one for runs you intend to publish.")
    return 2


def cmd_check(_args) -> int:
    """`refmasker check`: list each category's tools and versions, and the geNomad database.

    Returns 1 if a tool needed by a default category is missing.
    """
    from refmasker import database, external
    from refmasker.detectors import categories
    from refmasker.detectors.prophages import genomad_db, genomad_db_version

    print(f"refmasker {__version__}")
    missing_default = []
    for c in categories().values():
        # List both prophage callers: c.tools() needs the run options to pick one.
        tools = ["genomad", "PhiSpy.py"] if c.name == "prophage" else c.tools(None)
        for t in tools:
            v = external.tool_version(t)
            print(f"  {c.name:14} {t:12} {'NOT FOUND' if v == 'not found' else v}")
            # PhiSpy is optional (geNomad is the default caller), so it does not fail the check.
            if v == "not found" and c.default and t != "PhiSpy.py":
                missing_default.append(t)
    db = genomad_db(None)
    if db is None:
        print("  geNomad database: not set (use --genomad-db or REFMASKER_GENOMAD_DB)")
    else:
        print(f"  geNomad database: {db} (version {genomad_db_version(db)})")
        # Only read the recorded checksums; re-hashing 1.4 GB belongs in `download-db --verify`.
        mpath = database.manifest_path(db)
        if mpath.is_file():
            try:
                fingerprint = database.read_manifest(mpath).fingerprint or "not recorded"
            except (OSError, ValueError):
                fingerprint = f"unreadable ({mpath.name})"
            print(f"  database checksums: {fingerprint}")
        else:
            print(f"  database checksums: not recorded (no {mpath.name})")
    if missing_default:
        print(f"Missing tools for the default categories: {', '.join(missing_default)}")
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point of the `refmasker` command; returns the exit code.

    0 = success, 1 = no subcommand given or `check` found missing tools, 2 = error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 1
    setup_logging(getattr(args, "quiet", False))
    # Which function runs each subcommand.
    handlers = {
        "run": cmd_run,
        "compare": cmd_compare,
        "download-db": cmd_download_db,
        "check": cmd_check,
    }
    # Expected problems (bad input, unsafe or existing output folder, missing or failing tools;
    # ToolError is a RuntimeError) are logged as one clear message instead of a traceback.
    try:
        return handlers[args.command](args)
    except (ValueError, FileExistsError, FileNotFoundError, PermissionError, RuntimeError) as exc:
        logging.getLogger("refmasker").error("%s", exc)
        return 2
