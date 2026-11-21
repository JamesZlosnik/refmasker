"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from refmasker import __version__

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
        help="download the geNomad database (~1.4 GB)",
        description="Download the geNomad database used by the prophage category.",
    )
    d.add_argument("dest", help="folder to download into")

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


def cmd_download_db(args) -> int:
    """`refmasker download-db`: fetch the geNomad database with geNomad's own downloader."""
    from refmasker import external

    if external.missing_tools(["genomad"]):
        raise external.ToolError("genomad is not on PATH; install it from environment.yml first")
    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    external.run(["genomad", "download-database", str(dest)], dest, "genomad-download")
    print(f"geNomad database downloaded under {dest}. Point refmasker at the genomad_db folder")
    print(f"inside it, e.g.: export REFMASKER_GENOMAD_DB={dest / 'genomad_db'}")
    return 0


def cmd_check(_args) -> int:
    """`refmasker check`: list each category's tools and versions, and the geNomad database.

    Returns 1 if a tool needed by a default category is missing.
    """
    from refmasker import external
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
