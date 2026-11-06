from __future__ import annotations

"""
Frozen validated single-short-read caller.

Scientific definition inherited from the validated WGS_CONTA_2026 Ion Torrent
production branch:
  primary alignment, MAPQ >= 30, >= 80% query coverage, >= 98% identity,
  then coordinate-deduplicate passing single reads by (chromosome, left, right).

The helper functions cigar_stats() and parse_nm() are copied verbatim from the
validated production pipeline at patch-install time into
_single_legacy_analysis.py. This preserves the historical strict-count
definition exactly while making the installed SRAfeguard tree self-contained.
"""

import csv
import subprocess
from collections import Counter
from pathlib import Path
from ._single_legacy_analysis import cigar_stats, parse_nm


def _run(cmd, *, stdout=None, stderr=None):
    subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def human_map_single(index, reads, bam, log, threads, sort_threads):
    bam = Path(bam)
    log = Path(log)
    bam.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)

    with log.open("w") as lf:
        p1 = subprocess.Popen(
            ["minimap2", "-ax", "sr", "-t", str(threads), str(index), str(reads)],
            stdout=subprocess.PIPE,
            stderr=lf,
        )
        p2 = subprocess.Popen(
            ["samtools", "sort", "-@", str(sort_threads), "-o", str(bam)],
            stdin=p1.stdout,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        p1.stdout.close()
        _, err = p2.communicate()
        rc1 = p1.wait()

        if rc1 != 0 or p2.returncode != 0:
            raise RuntimeError(
                "Single-short human mapping failed "
                f"minimap2={rc1} samtools={p2.returncode}: "
                f"{(err or b'')[-1000:].decode(errors='ignore')}"
            )


def strict_human_single_from_bam(bam, evidence):
    evidence = Path(evidence)
    passing = []
    mapq30 = id95 = id98 = id99 = 0

    p = subprocess.Popen(
        ["samtools", "view", "-F", "0x904", str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    for line in p.stdout:
        f = line.rstrip().split("\t")
        qn = f[0]
        flag = int(f[1])
        chrom = f[2]
        pos = int(f[3])
        mapq = int(f[4])
        seq = f[9]

        if mapq < 30:
            continue

        mapq30 += 1
        q, a, r = cigar_stats(f[5])

        if not seq or not a:
            continue

        cov = 100 * q / len(seq)
        ident = 100 * (a - parse_nm(f)) / a

        if cov >= 80 and ident >= 95:
            id95 += 1
        if cov >= 80 and ident >= 98:
            id98 += 1
        if cov >= 80 and ident >= 99:
            id99 += 1

        if cov < 80 or ident < 98:
            continue

        passing.append((qn, chrom, pos, pos + r - 1, mapq, cov, ident, flag))

    _, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed: {err[-2000:]}")

    unique = {}
    for rec in passing:
        qn, chrom, left, right, *_ = rec
        key = (chrom, left, right)
        if key not in unique:
            unique[key] = rec

    evidence.mkdir(parents=True, exist_ok=True)

    qf = evidence / "strict_unique_qnames.txt"
    qf.write_text("".join(f"{v[0]}\n" for v in unique.values()))

    with (evidence / "strict_unique_reads.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["QNAME", "Chromosome", "Left", "Right", "MAPQ", "Coverage", "Identity"])
        for (_, _, _), (qn, chrom, left, right, mapq, cov, ident, flag) in unique.items():
            w.writerow([qn, chrom, left, right, mapq, f"{cov:.4f}", f"{ident:.4f}"])

    counts = Counter(k[0] for k in unique)

    with (evidence / "chromosome_counts.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Chromosome", "Strict_reads", "Fraction"])
        total = len(unique)
        for c, n in counts.most_common():
            w.writerow([c, n, f"{n / total:.6f}"])

    if counts:
        top_chr, top_n = counts.most_common(1)[0]
        top_frac = top_n / len(unique)
    else:
        top_chr, top_n, top_frac = "NA", 0, 0.0

    dist = "ZERO" if not counts else ("SINGLE_CHROMOSOME" if len(counts) == 1 else "MULTI_CHROMOSOME")

    return {
        "mapq30": mapq30,
        "id95": id95,
        "id98": id98,
        "id99": id99,
        "strict_reads": len(passing),
        "unique": len(unique),
        "qnames": qf,
        "nchr": len(counts),
        "distribution": dist,
        "top_chr": top_chr,
        "top_n": top_n,
        "top_frac": top_frac,
    }


def extract_strict_single(qnames, reads, evidence):
    evidence = Path(evidence)
    out = evidence / "strict_unique.fastq"
    qnames = Path(qnames)

    if (not qnames.exists()) or qnames.stat().st_size == 0:
        out.touch()
        print("+ [SRAfeguard fast-path] 0 strict QNAMEs; skipping seqkit FASTQ scan")
        return out

    _run(["seqkit", "grep", "-f", str(qnames), str(reads), "-o", str(out)])
    return out
