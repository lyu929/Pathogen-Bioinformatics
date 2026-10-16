"""Command-line interface: ``pathogenomics <command> ...``."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from pathogenomics import __version__


def _cmd_run(args: argparse.Namespace) -> int:
    from pathogenomics.config import load_config
    from pathogenomics.pipeline import STAGES, run_pipeline

    cfg = load_config(args.config)
    stages = args.only.split(",") if args.only else list(STAGES)
    run_pipeline(cfg, stages, write_report=not args.no_report)
    print(
        f"tables -> {cfg['output']['results']}\nreport -> {cfg['output']['report'] / 'REPORT.md'}"
    )
    return 0


def _cmd_trim(args: argparse.Namespace) -> int:
    from pathogenomics.qc import trim_fastq

    summary = trim_fastq(
        args.fastq, args.quality_cutoff, args.minimum_length, args.front_cutoff, args.output
    )
    print(json.dumps(summary.__dict__, indent=2))
    return 0


def _cmd_assembly(args: argparse.Namespace) -> int:
    from pathogenomics.assembly import assembly_table

    print(
        assembly_table(args.fasta, min_length=args.min_length, genome_size=args.genome_size).to_csv(
            sep="\t"
        ),
        end="",
    )
    return 0


def _cmd_fastq(args: argparse.Namespace) -> int:
    from pathogenomics.qc import fastq_stats

    print(json.dumps(fastq_stats(args.fastq).summary(), indent=2))
    return 0


def _cmd_pangenome(args: argparse.Namespace) -> int:
    from pathogenomics.pangenome import partition_summary, presence_matrix, read_roary

    print(partition_summary(presence_matrix(read_roary(args.csv))).to_csv(sep="\t"), end="")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pathogenomics", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="log progress to stderr")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("run", help="run the full pipeline and write results/ and reports/REPORT.md")
    p.add_argument("--config", type=Path, default=Path("config/course.toml"))
    p.add_argument(
        "--only", help="comma-separated stages: qc,assembly,annotation,pangenome,rnaseq,structure"
    )
    p.add_argument("--no-report", action="store_true", help="skip writing REPORT.md")
    p.set_defaults(func=_cmd_run)

    p = sub.add_parser(
        "trim", help="Cutadapt-compatible 3'/5' quality trimming of a single-end FASTQ"
    )
    p.add_argument("fastq", type=Path)
    p.add_argument("-q", "--quality-cutoff", type=int, default=20)
    p.add_argument("-m", "--minimum-length", type=int, default=22)
    p.add_argument("--front-cutoff", type=int, default=0)
    p.add_argument("-o", "--output", type=Path)
    p.set_defaults(func=_cmd_trim)

    p = sub.add_parser("fastq-stats", help="FastQC-style summary of a FASTQ")
    p.add_argument("fastq", type=Path)
    p.set_defaults(func=_cmd_fastq)

    p = sub.add_parser(
        "assembly-stats", help="N50/L50/GC/coverage for one or more FASTA assemblies"
    )
    p.add_argument("fasta", type=Path, nargs="+")
    p.add_argument("--min-length", type=int, default=0)
    p.add_argument("--genome-size", type=int, help="expected genome size for NG50")
    p.set_defaults(func=_cmd_assembly)

    p = sub.add_parser(
        "pangenome", help="core/shell/cloud partition of a Roary gene_presence_absence.csv"
    )
    p.add_argument("csv", type=Path)
    p.set_defaults(func=_cmd_pangenome)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
