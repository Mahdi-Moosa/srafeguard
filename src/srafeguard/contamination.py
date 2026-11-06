from pathlib import Path
import csv
import gzip
import json
import shutil
import tempfile

from .config import load_config
from . import paired_caller_frozen as paired_caller
from . import single_caller_frozen as single_caller


def _gzip_and_remove(path):
    path = Path(path)
    if not path.exists():
        return None
    gz = Path(str(path) + ".gz")
    with path.open("rb") as src, gzip.open(gz, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    path.unlink()
    return gz


def _count_data_rows(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with path.open() as fh:
        n = sum(1 for _ in fh)
    return max(0, n - 1)


def _count_nonempty_lines(path):
    path = Path(path)
    if not path.exists():
        return 0
    with path.open() as fh:
        return sum(1 for line in fh if line.strip())


def _chromosome_summary(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {
            "chromosomes_with_strict_evidence": None,
            "chromosome_counts": {},
        }

    counts = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames:
            chrom_candidates = ["Chromosome", "chromosome", "Contig", "contig"]
            count_candidates = [
                "Count", "count", "Strict_fragments", "strict_fragments",
                "Strict_reads", "strict_reads", "N", "n",
            ]
            chrom_col = next((x for x in chrom_candidates if x in reader.fieldnames), None)
            count_col = next((x for x in count_candidates if x in reader.fieldnames), None)
            if chrom_col and count_col:
                for row in reader:
                    try:
                        n = int(float(row[count_col]))
                    except Exception:
                        continue
                    counts[str(row[chrom_col])] = n

    return {
        "chromosomes_with_strict_evidence": (
            sum(v > 0 for v in counts.values()) if counts else None
        ),
        "chromosome_counts": counts,
    }


def run_contamination(
    fastq1,
    fastq2,
    output,
    sample,
    config_path=None,
    layout="PAIRED",
    platform="UNKNOWN",
):
    from .release_policy import validate_architecture
    validate_architecture(platform)
    cfg = load_config(config_path)
    fastq1 = Path(fastq1).resolve()
    fastq2 = Path(fastq2).resolve() if fastq2 else None
    output = Path(output).resolve()
    layout = str(layout).upper()

    if not fastq1.exists():
        raise FileNotFoundError(f"FASTQ1 not found: {fastq1}")

    if layout == "PAIRED":
        if fastq2 is None or not fastq2.exists():
            raise FileNotFoundError(f"FASTQ2 not found: {fastq2}")
    elif layout == "SINGLE":
        if fastq2 is not None:
            raise ValueError("SINGLE layout must not provide FASTQ2.")
    else:
        raise ValueError(f"Unsupported sequencing layout: {layout}")

    human_index = cfg["resources"]["human_minimap2_index"]
    if not human_index:
        raise RuntimeError("resources.human_minimap2_index is not configured")
    human_index = Path(human_index)
    if not human_index.exists():
        raise FileNotFoundError(f"Human minimap2 index not found: {human_index}")

    evidence = output / "evidence"
    logs = output / "logs"
    evidence.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    temp_root = cfg.get("runtime", {}).get("temp_dir")
    if temp_root:
        temp_root = Path(temp_root)
        temp_root.mkdir(parents=True, exist_ok=True)

    threads = int(
        cfg.get("runtime", {}).get(
            "mapper_threads", cfg.get("runtime", {}).get("threads", 8)
        )
    )
    sort_threads = int(cfg.get("runtime", {}).get("sort_threads", 2))
    keep_intermediates = bool(
        cfg.get("runtime", {}).get("keep_intermediates", False)
    )

    work = Path(
        tempfile.mkdtemp(
            prefix=f"hla_{sample}_",
            dir=str(temp_root) if temp_root else None,
        )
    )
    from .owned_temporary import register
    register(output,work,"mapping_workspace")
    bam = work / f"{sample}.human.bam"
    mmap_log = logs / "human_minimap2.log"

    try:
        if layout == "PAIRED":
            paired_caller.human_map(
                human_index, fastq1, fastq2, bam, mmap_log, threads, sort_threads
            )
            paired_caller.strict_human_from_bam(bam, evidence)
            qnames = evidence / "strict_unique_qnames.txt"
            paired_caller.extract_strict(qnames, fastq1, fastq2, evidence)

            strict_r1 = _gzip_and_remove(evidence / "strict_unique_R1.fastq")
            strict_r2 = _gzip_and_remove(evidence / "strict_unique_R2.fastq")
            strict_table = evidence / "strict_unique_fragments.tsv"
            evidence_unit = "physical_pair"
        else:
            single_caller.human_map_single(
                human_index, fastq1, bam, mmap_log, threads, sort_threads
            )
            single_caller.strict_human_single_from_bam(bam, evidence)
            qnames = evidence / "strict_unique_qnames.txt"
            single_caller.extract_strict_single(qnames, fastq1, evidence)

            strict_r1 = _gzip_and_remove(evidence / "strict_unique.fastq")
            strict_r2 = None
            strict_table = evidence / "strict_unique_reads.tsv"
            evidence_unit = "single_read"

        strict_fragments = _count_data_rows(strict_table)
        strict_qnames = _count_nonempty_lines(qnames)
        chrom = _chromosome_summary(evidence / "chromosome_counts.tsv")

        summary = {
            "sample": sample,
            "layout": layout,
            "platform": platform or "UNKNOWN",
            "evidence_unit": evidence_unit,
            "fastq1": str(fastq1),
            "fastq2": str(fastq2) if fastq2 else None,
            "strict_human_fragments": strict_fragments,
            "strict_unique_qnames": strict_qnames,
            "chromosomes_with_strict_evidence": chrom["chromosomes_with_strict_evidence"],
            "strict_R1": str(strict_r1) if strict_r1 else None,
            "strict_R2": str(strict_r2) if strict_r2 else None,
            "human_mapping_log": str(mmap_log),
        }

        with (output / "contamination_summary.json").open("w") as fh:
            json.dump(summary, fh, indent=2)

        with (output / "contamination_summary.tsv").open("w") as fh:
            fields = [
                "sample", "layout", "platform", "evidence_unit",
                "strict_human_fragments", "strict_unique_qnames",
                "chromosomes_with_strict_evidence", "strict_R1", "strict_R2",
            ]
            fh.write("\t".join(fields) + "\n")
            fh.write(
                "\t".join(
                    "" if summary[x] is None else str(summary[x]) for x in fields
                ) + "\n"
            )

        return summary
    finally:
        if not keep_intermediates:
            shutil.rmtree(work, ignore_errors=True)
