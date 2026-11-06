from pathlib import Path
import re
import argparse
import sys

from .runner import run_audit
from .runtime_config import resolve_config_path
from .resources import print_resource_check
from .batch import run_batch
from .reanalysis import run_reanalysis


def build_parser():
    parser=argparse.ArgumentParser(
        prog="srafeguard",
        description=(
            "Audit recoverable human genetic information "
            "in nominally microbial sequencing datasets."
        ),
    )

    sub=parser.add_subparsers(
        dest="command",
        required=True,
    )

    resources=sub.add_parser(
        "resources",
        help="Manage reference resources.",
    )

    rsub=resources.add_subparsers(
        dest="resources_command",
        required=True,
    )

    check=rsub.add_parser(
        "check",
        help="Check reference resources and external tools.",
    )

    check.add_argument("--config")

    audit=sub.add_parser(
        "audit",
        help="Audit one sequencing library for recoverable human genomic information.",
    )

    group=audit.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "--sra",
        help="Local SRA file.",
    )

    group.add_argument(
        "--fastq1",
        help="Read-1 FASTQ/FASTQ.gz.",
    )

    audit.add_argument(
        "--fastq2",
        help="Read-2 FASTQ/FASTQ.gz. Omit for single-end input.",
    )

    audit.add_argument(
        "--layout",
        choices=["auto", "paired", "single"],
        default="auto",
        help=(
            "Sequencing layout. Default: auto. For SRA input, layout is "
            "detected from SRA READ_TYPE metadata; for FASTQ input, "
            "presence/absence of --fastq2 is used."
        ),
    )

    audit.add_argument(
        "--platform",
        default="auto",
        help=(
            "Platform override. Default: auto. For SRA input, SRA metadata "
            "is used when available. Caller selection is architecture-first."
        ),
    )

    audit.add_argument(
        "--sample",
        help=(
            "Sample identifier. For SRA input this is "
            "normally inferred from the filename."
        ),
    )

    audit.add_argument(
        "--reject-source-sra",
        help=(
            "For externally filtered FASTQs, use this SRA "
            "to apply the WGS_CONTA_2026 biological "
            "READ_FILTER policy before analysis."
        ),
    )

    audit.add_argument(
        "--output",
    )

    audit.add_argument(
        "--config",
    )

    audit.add_argument(
        "--keep-evidence",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Retain strict-human BAM/FASTQs/QNAME evidence after analysis; --no-keep-evidence forces cleanup.",
    )


    batch=sub.add_parser(
        "batch",
        help="Run a resume-safe TSV manifest of audits.",
    )
    batch.add_argument("--manifest", required=True)
    batch.add_argument("--output-root", required=True)
    batch.add_argument("--config")
    batch.add_argument("--keep-evidence", action=argparse.BooleanOptionalAction, default=None)
    batch.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    batch.add_argument("--fail-fast", action="store_true")

    rea=sub.add_parser(
        "reanalyze",
        help="Re-run a completed result using its saved effective configuration as the baseline.",
    )
    rea.add_argument("--results", required=True)
    rea.add_argument("--output")
    rea.add_argument("--config", help="Current configuration overlay on inherited settings.")
    rea.add_argument("--min-mapq", type=int)
    rea.add_argument("--min-bq", type=int)
    rea.add_argument("--set", action="append", default=[], metavar="DOTTED.KEY=VALUE")
    rea.add_argument("--keep-evidence", action=argparse.BooleanOptionalAction, default=None)

    return parser


def _infer_sample_name(args):
    """Infer a stable sample name for default output paths."""
    if getattr(args, "sample", None):
        return str(args.sample)

    sra = getattr(args, "sra", None)
    if sra:
        name = Path(sra).name
        if name.lower().endswith(".sra"):
            name = name[:-4]
        return name or "sample"

    fq1 = getattr(args, "fastq1", None)
    if fq1:
        name = Path(fq1).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.lower().endswith(suffix):
                name = name[:-len(suffix)]
                break
        name = re.sub(r"(?i)(?:[_\.-]R?1)$", "", name)
        return name or "sample"

    return "sample"


def _default_output_for_args(args):
    """Default per-sample output directory."""
    sample = _infer_sample_name(args)
    return str(Path.cwd() / "srafeguard_results" / sample)


def main():
    import os
    os.umask(0o077)
    parser=build_parser()

    args=parser.parse_args()

    # SRAFEGUARD_AUTO_CONFIG_DISCOVERY
    if hasattr(args, "config"):
        try:
            args.config = resolve_config_path(args.config)
        except Exception as exc:
            raise SystemExit(f"SRAfeguard error: {exc}") from None
    if (
        args.command=="resources"
        and args.resources_command=="check"
    ):
        return print_resource_check(
            args.config
        )

    if args.command=="audit":
        if args.fastq1 and args.layout == "paired" and not args.fastq2:
            parser.error("--fastq2 is required with --layout paired")
        if args.fastq1 and args.layout == "single" and args.fastq2:
            parser.error("--fastq2 cannot be used with --layout single")

        result=run_audit(
            sra=args.sra,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            sample=args.sample,
            reject_source_sra=args.reject_source_sra,
            layout=args.layout,
            platform=args.platform,
            output=args.output or _default_output_for_args(args),
            config_path=args.config,
            keep_evidence=args.keep_evidence,
        )

        print()
        print("="*80)
        print("SRAFEGUARD")
        print("="*80)
        print(
            f"Sample:                  "
            f"{result['sample']}"
        )
        print(
            f"Strict human fragments:  "
            f"{result['strict_human_fragments']:,}"
        )
        print(
            f"Strict unique QNAMEs:     "
            f"{result['strict_unique_qnames']:,}"
        )
        print(
            f"Chromosomes represented:  "
            f"{result['chromosomes_with_strict_evidence']}"
        )

        return 0


    if args.command=="batch":
        run_batch(
            manifest=args.manifest,
            output_root=args.output_root,
            config_path=args.config,
            keep_evidence=args.keep_evidence,
            resume=args.resume,
            fail_fast=args.fail_fast,
        )
        return 0

    if args.command=="reanalyze":
        run_reanalysis(
            results=args.results,
            output=args.output,
            config_path=args.config,
            min_mapq=args.min_mapq,
            min_bq=args.min_bq,
            sets=args.set,
            keep_evidence=args.keep_evidence,
        )
        return 0

    parser.error("Unsupported command")


if __name__=="__main__":
    sys.exit(main())
