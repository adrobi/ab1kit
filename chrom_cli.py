"""ab1kit entry point: the GUI by default, plus ``verify-clones`` / ``qc`` / ``export``.

The subcommands are a thin shell over ``chrom_core`` and never import Qt.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

COMMANDS = {"verify-clones", "qc", "export"}


class CliError(Exception):
    """Bad input: reported on stderr, exit code 2."""


def _console_encoding(code_page: int) -> str:
    """Python codec for a console code page; utf-8 when unknown (0 = query failed)."""
    import codecs

    if code_page:
        try:
            codecs.lookup(f"cp{code_page}")
            return f"cp{code_page}"
        except LookupError:
            pass
    return "utf-8"


def _attach_console() -> None:
    """Frozen windowed build: reuse the parent console, writing in its own code page."""
    opened = []
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        if not kernel32.AttachConsole(-1):
            return
        encoding = _console_encoding(kernel32.GetConsoleOutputCP())  # never change it
        wanted = [name for name in ("stdout", "stderr") if getattr(sys, name) is None]
        for name in wanted:
            opened.append((name, open("CONOUT$", "w", encoding=encoding, errors="replace")))
    except (AttributeError, OSError, LookupError):
        for _name, stream in opened:  # all or nothing: keep both streams consistent
            try:
                stream.close()
            except (AttributeError, OSError):
                pass
        return
    for name, stream in opened:
        setattr(sys, name, stream)


def _prepare_streams() -> None:
    """Print Cyrillic safely on any console."""
    attached = sys.stdout is None or sys.stderr is None
    if attached:
        _attach_console()
    for stream in (sys.stdout, sys.stderr):
        # an attached console stream keeps the console's own code page
        if stream is not None and hasattr(stream, "reconfigure") and not (
                attached and getattr(stream, "name", "") == "CONOUT$"):
            stream.reconfigure(encoding="utf-8", errors="replace")


def _parse_region(text: str) -> tuple[int, int]:
    try:
        a, b = (int(x) for x in text.split("-"))
    except ValueError:
        raise CliError(f"invalid --region {text!r}: expected A-B (1-based, inclusive)") from None
    if a < 1 or b < 1:
        raise CliError(f"invalid --region {text!r}: positions start at 1")
    return (a - 1, b)


def _reads_in(folder: str) -> list[Path]:
    path = Path(folder)
    if not path.is_dir():
        raise CliError(f"not a directory: {folder}")
    return sorted(p for p in path.iterdir() if p.suffix.lower() == ".ab1")


def _verify_clones(args: argparse.Namespace) -> int:
    from chrom_core.clone_report import (
        clone_report_xlsx, clones_csv, differences_csv, read_issue)
    from chrom_core.clones import VerifyParams, check_region, group_reads, verify_all
    from chrom_core.reference import load_construct

    try:
        construct = load_construct(args.construct)
    except (ValueError, OSError) as exc:
        raise CliError(f"cannot load construct: {exc}") from exc
    if args.linear:
        construct = construct.with_circular(False)
    out = Path(args.output)
    if out.suffix.lower() not in (".xlsx", ".csv"):
        raise CliError("report file must end in .xlsx or .csv")
    region = _parse_region(args.region) if args.region else None
    try:
        check_region(region, len(construct.seq), construct.circular)
        paths = _reads_in(args.reads_dir)
        if not paths:
            raise CliError(f"no .ab1 files in {args.reads_dir}")
        grouping = group_reads(paths, pattern=args.pattern)
    except ValueError as exc:
        raise CliError(str(exc)) from exc
    for path in grouping.unassigned:
        print(f"ab1kit: skipped (clone name not recognised): {path.name}", file=sys.stderr)
    if not grouping.reads:
        raise CliError("no reads could be assigned to a clone")
    params = VerifyParams(min_q=args.min_q, het_ratio=args.het_ratio)
    reports = verify_all(construct, region, grouping, params)
    for report in reports:
        for read in report.reads:
            if read_issue(read):
                print(f"ab1kit: {read.ref.path.name}: {read_issue(read)}", file=sys.stderr)

    if out.suffix.lower() == ".xlsx":
        clone_report_xlsx(out, reports, construct.name)
    else:
        out.write_text(clones_csv(reports), encoding="utf-8")
        out.with_name(f"{out.stem}_differences.csv").write_text(
            differences_csv(reports), encoding="utf-8")

    width = max(len(r.clone) for r in reports)
    for r in reports:
        print(f"{r.clone:<{width}}  {r.verdict:<10}  coverage {r.region_coverage_pct:5.1f}%  "
              f"differences {len(r.error_differences)}")
    return 0 if any(r.verdict == "ok" for r in reports) else 1


def _qc(args: argparse.Namespace) -> int:
    from chrom_core import batch

    _reads_in(args.reads_dir)  # validates the folder
    text = batch.summary_csv(batch.scan_folder(args.reads_dir))
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


def _export(args: argparse.Namespace) -> int:
    from chrom_core import batch

    if not args.trim:
        raise CliError("export currently needs --trim")
    reports = batch.scan_folder(args.reads_dir)
    if not reports:
        raise CliError(f"no .ab1 files in {args.reads_dir}")
    results = batch.batch_export(reports, args.output, fmt="fasta" if args.fasta else "fastq")
    failed = 0
    for source, _dest, error in results:
        if error:
            failed += 1
            print(f"ab1kit: {source.name}: {error}", file=sys.stderr)
    print(f"exported {len(results) - failed} of {len(results)} reads to {args.output}")
    return 1 if failed else 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ab1kit")
    sub = parser.add_subparsers(dest="command", required=True)

    v = sub.add_parser("verify-clones", help="verify clones against a construct")
    v.add_argument("--construct", required=True, help="FASTA, GenBank or SnapGene file")
    v.add_argument("--region", help="A-B, 1-based inclusive (default: whole construct)")
    v.add_argument("--linear", action="store_true", help="treat the construct as linear")
    v.add_argument("--pattern", help="regex with named groups 'clone' and 'primer'")
    v.add_argument("--min-q", type=int, default=20)
    v.add_argument("--het-ratio", type=float, default=0.40)
    v.add_argument("reads_dir", metavar="READS_DIR")
    v.add_argument("-o", "--output", required=True, help="report: .xlsx or .csv")
    v.set_defaults(func=_verify_clones)

    q = sub.add_parser("qc", help="quality summary of every .ab1 in a folder")
    q.add_argument("reads_dir", metavar="READS_DIR")
    q.add_argument("-o", "--output", help="CSV file (default: stdout)")
    q.set_defaults(func=_qc)

    e = sub.add_parser("export", help="export trimmed reads")
    e.add_argument("reads_dir", metavar="READS_DIR")
    e.add_argument("--trim", action="store_true")
    fmt = e.add_mutually_exclusive_group()
    fmt.add_argument("--fasta", action="store_true")
    fmt.add_argument("--fastq", action="store_true", help="(default)")
    e.add_argument("-o", "--output", required=True, metavar="DIR")
    e.set_defaults(func=_export)
    return parser


def _run_command(argv: list[str]) -> int:
    _prepare_streams()
    try:
        args = _build_parser().parse_args(argv)
    except SystemExit as exc:  # argparse usage errors and --help
        return exc.code if isinstance(exc.code, int) else 2
    try:
        return args.func(args)
    except CliError as exc:
        print(f"ab1kit: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ab1kit: {exc}", file=sys.stderr)
        return 2


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in COMMANDS:
        return _run_command(args)
    from ui.main_window import run

    return run(args)


if __name__ == "__main__":
    sys.exit(main())
