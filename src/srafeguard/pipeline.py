from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .alignment_evidence import (
    capture_original_human_bam,
    finalize_strict_human_bam,
)
from .audit import run_audit as _run_strict_audit
from .clinvar import run_clinvar
from .single_clinvar import run_clinvar_single
from .population import run_population
from .runtime_config import load_effective_config, resolve_config_path


def _infer_sample(
    *,
    sample: str | None,
    sra: str | None,
    fastq1: str | None,
) -> str:
    if sample:
        return str(sample)

    if sra:
        name = Path(sra).name
        return name[:-4] if name.lower().endswith(".sra") else name

    if fastq1:
        name = Path(fastq1).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.lower().endswith(suffix):
                name = name[:-len(suffix)]
                break
        for suffix in ("_R1", "_1", ".R1", ".1", "-R1", "-1"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        return name

    return "sample"


def _input_label(
    *,
    sra: str | None,
    fastq1: str | None,
    fastq2: str | None,
) -> str:
    if sra:
        return Path(sra).name

    if fastq1 and fastq2:
        return f"{Path(fastq1).name} + {Path(fastq2).name}"

    if fastq1:
        return Path(fastq1).name

    return "unknown"


def _validate_inputs(
    *,
    sra: str | None,
    fastq1: str | None,
    fastq2: str | None,
    reject_source_sra: str | None,
):
    for label, value in (
        ("SRA input", sra),
        ("FASTQ R1", fastq1),
        ("FASTQ R2", fastq2),
        ("reject-source SRA", reject_source_sra),
    ):
        if value and not Path(value).expanduser().is_file():
            raise FileNotFoundError(f"{label} not found: {value}")


def _tool(cfg: dict, key: str, fallback: str) -> str:
    tools = cfg.get("tools", {}) or {}
    value = tools.get(key)

    if value:
        value = str(value)

        if "/" in value:
            if not Path(value).exists():
                raise FileNotFoundError(
                    f"Configured tool not found: {value}"
                )
            return value

        found = shutil.which(value)
        if found:
            return found

    found = shutil.which(fallback)
    if not found:
        raise FileNotFoundError(
            f"Required tool '{fallback}' not found in PATH."
        )

    return found


def _runtime_root(cfg: dict) -> Path:
    runtime = cfg.get("runtime", {}) or {}
    value = runtime.get("tempdir") or runtime.get("temp_dir")

    if value:
        root = Path(value)
    else:
        root = Path(tempfile.gettempdir()) / "srafeguard"

    root.mkdir(parents=True, exist_ok=True)
    return root


def _read_one_row(path: Path) -> dict:
    if not path.is_file():
        return {}

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}

    return df.iloc[0].to_dict()


def _pick(d: dict, *keys, default=None):
    for key in keys:
        if key in d and pd.notna(d[key]):
            return d[key]
    return default


def _write_privacy_summary(
    *,
    output: Path,
    sample: str,
    input_label: str,
    population: dict,
    clinvar: dict,
):
    contamination = _read_one_row(
        output / "contamination_summary.tsv"
    )

    row = {
        "Input": input_label,
        "Sample": sample,
        "Layout": _pick(contamination, "layout", "Layout", default=""),
        "Platform": _pick(contamination, "platform", "Platform", default=""),
        "Evidence_unit": _pick(contamination, "evidence_unit", "Evidence_unit", default=""),
        "ClinVar_evaluation_status": clinvar.get(
            "ClinVar_evaluation_status", "COMPLETE"
        ),
        "ClinVar_validation_mode": clinvar.get(
            "ClinVar_validation_mode",
            "paired_physical_pair" if str(
                _pick(contamination, "layout", "Layout", default="")
            ).upper() == "PAIRED" else "",
        ),
        "ClinVar_validated_P_LP_supporting_reads": clinvar.get(
            "ClinVar_validated_P_LP_supporting_reads", ""
        ),
        "Strict_human_fragments": _pick(
            contamination,
            "strict_human_fragments",
            "Strict_human_fragments",
            default="",
        ),
        "Strict_unique_QNAMEs": _pick(
            contamination,
            "strict_unique_qnames",
            "Strict_unique_QNAMEs",
            default="",
        ),
        "Human_chromosomes_represented": _pick(
            contamination,
            "chromosomes_with_strict_evidence",
            "Chromosomes_represented",
            default="",
        ),
        "Population_sites_used": population.get("Sites_used", 0),
        "Population_distinct_1Mb_blocks": population.get(
            "Distinct_1Mb_blocks", 0
        ),
        "Population_best_affinity": population.get(
            "Best_population", ""
        ),
        "Population_second_affinity": population.get(
            "Second_population", ""
        ),
        "Population_best_bootstrap": population.get(
            "Best_bootstrap_support", 0.0
        ),
        "Population_affinity_recoverable": population.get(
            "Population_affinity_recoverable", False
        ),
        "Raw_ClinVar_ALT_sites_observed": clinvar.get(
            "Raw_ClinVar_ALT_sites_observed", 0
        ),
        "Raw_ClinVar_candidate_P_LP_ALT_sites": clinvar.get(
            "Raw_ClinVar_candidate_P_LP_ALT_sites", 0
        ),
        "Hotspot_ClinVar_ALT_sites_observed": clinvar.get(
            "Hotspot_ClinVar_ALT_sites_observed", 0
        ),
        "Hotspot_candidate_P_LP_ALT_sites": clinvar.get(
            "Hotspot_candidate_P_LP_ALT_sites", 0
        ),
        "ClinVar_sites_with_HQ_base_observation": clinvar.get(
            "ClinVar_sites_with_HQ_base_observation", 0
        ),
        "ClinVar_ALT_sites_observed": clinvar.get(
            "ClinVar_ALT_sites_observed", 0
        ),
        "ClinVar_candidate_P_LP_ALT_sites": clinvar.get(
            "ClinVar_candidate_P_LP_ALT_sites", 0
        ),
        "ClinVar_validated_P_LP_ALT_sites": clinvar.get(
            "ClinVar_validated_P_LP_ALT_sites", 0
        ),
        "ClinVar_validated_P_LP_supporting_mates": clinvar.get(
            "ClinVar_validated_P_LP_supporting_mates", 0
        ),
        "ClinVar_exposure_detected": clinvar.get(
            "ClinVar_exposure_detected", False
        ),
        "P_LP_candidate_exposure_detected": clinvar.get(
            "P_LP_candidate_exposure_detected", False
        ),
        "P_LP_validated_exposure_detected": clinvar.get(
            "P_LP_validated_exposure_detected", False
        ),
    }

    pd.DataFrame([row]).to_csv(
        output / "privacy_summary.tsv",
        sep="\t",
        index=False,
    )


def _cleanup_read_evidence(output: Path):
    evidence = output / "evidence"

    if not evidence.exists():
        return

    for name in (
        "strict_human.bam",
        "strict_human.bam.bai",
        "strict_unique_R1.fastq.gz",
        "strict_unique_R2.fastq.gz",
        "strict_unique_R1.fastq",
        "strict_unique_R2.fastq",
        "strict_unique.fastq.gz",
        "strict_unique.fastq",
        "strict_unique_qnames.txt",
    ):
        (evidence / name).unlink(missing_ok=True)

    try:
        evidence.rmdir()
    except OSError:
        pass


def run_audit(
    *,
    sra=None,
    fastq1=None,
    fastq2=None,
    sample=None,
    reject_source_sra=None,
    layout=None,
    platform=None,
    output=None,
    config_path=None,
    keep_evidence=None,
    **kwargs,
):
    """Complete public SRAfeguard audit."""
    captured_bam = None

    try:
        _validate_inputs(
            sra=sra,
            fastq1=fastq1,
            fastq2=fastq2,
            reject_source_sra=reject_source_sra,
        )

        resolved_config = resolve_config_path(config_path)
        cfg, _ = load_effective_config(resolved_config)

        resolved_sample = _infer_sample(
            sample=sample,
            sra=sra,
            fastq1=fastq1,
        )

        input_label = _input_label(
            sra=sra,
            fastq1=fastq1,
            fastq2=fastq2,
        )

        if keep_evidence is None:
            keep_evidence = bool(
                (cfg.get("outputs", {}) or {}).get(
                    "keep_evidence", False
                )
            )

        print()
        print("=" * 80)
        print("SRAFEGUARD — INPUT")
        print("=" * 80)
        print(f"Input:   {input_label}")
        print(f"Sample:  {resolved_sample}")
        print(
            "Evidence retention: "
            + ("KEEP" if keep_evidence else "CLEAN AFTER ANALYSIS")
        )
        print()

        output_path = Path(output).resolve()
        evidence = output_path / "evidence"
        evidence.mkdir(parents=True, exist_ok=True)

        runtime_root = _runtime_root(cfg)
        capture_dir = runtime_root / "captured_original_bams"
        capture_dir.mkdir(parents=True, exist_ok=True)

        captured_bam = capture_dir / (
            f"{resolved_sample}.original_human_mapping.bam"
        )

        from .owned_temporary import register
        register(output_path,captured_bam,"captured_mapping")
        with capture_original_human_bam(captured_bam):
            result = _run_strict_audit(
                sra=sra,
                fastq1=fastq1,
                fastq2=fastq2,
                sample=sample,
                reject_source_sra=reject_source_sra,
                layout=layout,
                platform=platform,
                output=output,
                config_path=resolved_config,
                **kwargs,
            )

        samtools = _tool(cfg, "samtools", "samtools")
        sort_threads = int(
            (cfg.get("runtime", {}) or {}).get(
                "sort_threads", 2
            )
        )

        strict_bam = evidence / "strict_human.bam"
        qnames = evidence / "strict_unique_qnames.txt"

        finalize_strict_human_bam(
            captured_bam=captured_bam,
            qnames_file=qnames,
            output_bam=strict_bam,
            samtools=samtools,
            sort_threads=sort_threads,
        )

        captured_bam.unlink(missing_ok=True)

        population = run_population(
            output=output_path,
            sample=resolved_sample,
            config_path=resolved_config,
        )

        resolved_layout = str(result.get("layout", "PAIRED")).upper()

        if resolved_layout == "PAIRED":
            clinvar = run_clinvar(
                output=output_path,
                sample=resolved_sample,
                config_path=resolved_config,
            )
            clinvar["ClinVar_evaluation_status"] = "COMPLETE"
        else:
            clinvar = run_clinvar_single(
                output=output_path,
                sample=resolved_sample,
                config_path=resolved_config,
            )

        _write_privacy_summary(
            output=output_path,
            sample=resolved_sample,
            input_label=input_label,
            population=population,
            clinvar=clinvar,
        )

        if not keep_evidence:
            _cleanup_read_evidence(output_path)

        return result

    except KeyboardInterrupt:
        raise
    except SystemExit:
        raise
    except subprocess.CalledProcessError as exc:
        cmd = (
            " ".join(map(str, exc.cmd))
            if exc.cmd
            else "external command"
        )
        raise SystemExit(
            f"SRAfeguard error: command failed with exit code "
            f"{exc.returncode}: {cmd}"
        ) from None
    except Exception as exc:
        raise SystemExit(f"SRAfeguard error: {exc}") from None
    finally:
        if captured_bam is not None:
            Path(captured_bam).unlink(missing_ok=True)
