from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from .population import base_at_reference_positions
from .runtime_config import load_effective_config


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _cfg_get(d: dict, *keys, default=None):
    cur: Any = d
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _tool(cfg: dict, key: str, fallback: str) -> str:
    value = _cfg_get(cfg, "tools", key, default=None)

    if value:
        value = str(value)
        if os.path.sep in value:
            if not Path(value).exists():
                raise FileNotFoundError(f"Configured tool not found: {value}")
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


def _settings(cfg: dict) -> dict:
    c = cfg.get("clinvar", {}) or {}

    validation = c.get("validation", {}) or {}

    return {
        "enabled": bool(c.get("enabled", True)),
        "min_mapq": int(c.get("min_mapq", c.get("mapq", 30))),
        "min_bq": int(c.get("min_bq", c.get("bq", 30))),
        "biallelic_snv_only": bool(
            c.get("biallelic_snv_only", c.get("biallelic_snv", True))
        ),
        "apply_hotspot_mask": bool(
            c.get(
                "apply_hotspot_exclusion_mask",
                c.get("apply_hotspot_mask", True),
            )
        ),
        "chromosome_scope": str(c.get("chromosome_scope", "all")).lower(),
        "validate_plp": bool(validation.get("enabled", True)),
        "validation_threads": int(validation.get("threads", 4)),
        "validation_secondary_n": int(validation.get("secondary_n", 50)),
        "validation_secondary_p": float(validation.get("secondary_p", 0.50)),
        "validation_min_remap_mapq": int(
            validation.get("min_remap_mapq", 30)
        ),
        "validation_distinct_locus_window_bp": int(
            validation.get("distinct_locus_window_bp", 500)
        ),
        "validation_min_passing_mates": int(
            validation.get("min_passing_mates", 1)
        ),
    }


def _reference_end(start_1based: int, cigar: str) -> int:
    """Inclusive reference end of an alignment."""
    rpos = start_1based
    consumed = 0

    for ns, op in CIGAR_RE.findall(cigar):
        n = int(ns)
        if op in {"M", "=", "X", "D", "N"}:
            consumed += n

    if consumed == 0:
        return start_1based

    return start_1based + consumed - 1


def _merge_bed(intervals):
    """Merge 0-based half-open BED intervals."""
    by_chrom = defaultdict(list)

    for chrom, start, end in intervals:
        start = int(start)
        end = int(end)
        if end <= start:
            continue
        by_chrom[str(chrom)].append((start, end))

    out = []

    for chrom, vals in by_chrom.items():
        vals.sort()
        cur_s, cur_e = vals[0]

        for s, e in vals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e

        out.append((chrom, cur_s, cur_e))

    return out


def _read_bed(path: Path):
    rows = []

    if not path.is_file():
        return rows

    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue

            rows.append((f[0], int(f[1]), int(f[2])))

    return _merge_bed(rows)


def _subtract_bed(intervals, masks):
    """Subtract mask BED intervals from BED intervals."""
    masks_by_chrom = defaultdict(list)
    for chrom, s, e in _merge_bed(masks):
        masks_by_chrom[chrom].append((s, e))

    result = []

    for chrom, start, end in _merge_bed(intervals):
        segments = [(start, end)]

        for ms, me in masks_by_chrom.get(chrom, []):
            if me <= start:
                continue
            if ms >= end:
                break

            new_segments = []

            for s, e in segments:
                if me <= s or ms >= e:
                    new_segments.append((s, e))
                    continue

                if ms > s:
                    new_segments.append((s, min(ms, e)))
                if me < e:
                    new_segments.append((max(me, s), e))

            segments = new_segments
            if not segments:
                break

        for s, e in segments:
            if e > s:
                result.append((chrom, s, e))

    return result


def _strict_fragment_territory(
    *,
    bam: Path,
    samtools: str,
    min_mapq: int,
):
    """Derive physical strict-fragment spans from the retained BAM.

    Paired alignments sharing a QNAME are collapsed to the outer genomic span.
    This mirrors the privacy question used in the study: what ClinVar loci lie
    within territory physically represented by a strict human fragment?
    """
    qspan = {}

    p = subprocess.Popen(
        [samtools, "view", "-F", "2308", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert p.stdout is not None

    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue

        qname = f[0]
        chrom = f[2]
        if chrom == "*":
            continue

        try:
            mapq = int(f[4])
            start = int(f[3])
        except ValueError:
            continue

        if mapq < min_mapq:
            continue

        cigar = f[5]
        if cigar == "*":
            continue

        end = _reference_end(start, cigar)
        key = (qname, chrom)

        if key not in qspan:
            qspan[key] = [start, end]
        else:
            qspan[key][0] = min(qspan[key][0], start)
            qspan[key][1] = max(qspan[key][1], end)

    rc = p.wait()
    if rc != 0:
        raise RuntimeError(
            f"samtools view failed while deriving ClinVar territory: {rc}"
        )

    # Convert inclusive 1-based span to 0-based half-open BED.
    intervals = [
        (chrom, start - 1, end)
        for (_, chrom), (start, end) in qspan.items()
    ]

    return _merge_bed(intervals)


def _write_bed(intervals, path: Path):
    with path.open("w") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def _parse_info(info_text: str) -> dict:
    info = {}

    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True

    return info


def _is_conservative_plp(clnsig: str) -> bool:
    """Conservative P/LP annotation flag.

    Conflicting/uncertain/benign-containing records are not counted as
    canonical P/LP exposure even if a pathogenic token also appears.
    """
    s = (clnsig or "").lower()

    disqualifiers = (
        "conflict",
        "uncertain",
        "benign",
        "drug_response",
        "risk_factor",
        "protective",
    )

    if any(x in s for x in disqualifiers):
        return False

    return (
        "pathogenic" in s
        or "likely_pathogenic" in s
        or "likely pathogenic" in s
    )


def _query_clinvar(
    *,
    vcf: Path,
    territory_bed: Path,
    bcftools: str,
    biallelic_snv_only: bool,
):
    cmd = [bcftools, "view", "-R", str(territory_bed)]

    if biallelic_snv_only:
        cmd += ["-m2", "-M2", "-v", "snps"]

    cmd += [str(vcf)]

    p = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )

    rows = []

    for line in p.stdout.splitlines():
        if not line or line.startswith("#"):
            continue

        f = line.split("\t")
        if len(f) < 8:
            continue

        chrom, pos, vid, ref, alt, qual, flt, info_text = f[:8]

        if biallelic_snv_only:
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref.upper() not in "ACGT" or alt.upper() not in "ACGT":
                continue

        info = _parse_info(info_text)
        clnsig = str(info.get("CLNSIG", ""))

        rows.append(
            {
                "Chromosome": chrom,
                "Position": int(pos),
                "ID": vid,
                "REF": ref.upper(),
                "ALT": alt.upper(),
                "CLNSIG": clnsig,
                "CLNREVSTAT": str(info.get("CLNREVSTAT", "")),
                "CLNDN": str(info.get("CLNDN", "")),
                "GENEINFO": str(info.get("GENEINFO", "")),
                "P_LP_annotation": _is_conservative_plp(clnsig),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "Chromosome", "Position", "ID", "REF", "ALT",
                "CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO",
                "P_LP_annotation",
            ]
        )

    return (
        pd.DataFrame(rows)
        .sort_values(["Chromosome", "Position", "REF", "ALT"])
        .drop_duplicates(
            ["Chromosome", "Position", "REF", "ALT"],
            keep="first",
        )
    )


def _observe_candidate_sites(
    *,
    bam: Path,
    candidates: pd.DataFrame,
    samtools: str,
    min_mapq: int,
    min_bq: int,
):
    if candidates.empty:
        return candidates.assign(
            High_quality_base_observations=pd.Series(dtype=int),
            REF_observations=pd.Series(dtype=int),
            ALT_observations=pd.Series(dtype=int),
            Other_observations=pd.Series(dtype=int),
            Mean_BQ=pd.Series(dtype=float),
            Mean_MAPQ=pd.Series(dtype=float),
            ALT_actually_observed=pd.Series(dtype=bool),
            P_LP_ALT_observed=pd.Series(dtype=bool),
        )

    by_chrom = {}
    lookup = {}

    for row in candidates.itertuples(index=False):
        chrom = str(row.Chromosome)
        pos = int(row.Position)

        by_chrom.setdefault(chrom, []).append(pos)
        lookup[(chrom, pos, row.REF, row.ALT)] = row

    for chrom in by_chrom:
        by_chrom[chrom] = sorted(set(by_chrom[chrom]))

    counts = defaultdict(lambda: defaultdict(int))
    bq_sum = defaultdict(int)
    mapq_sum = defaultdict(int)

    p = subprocess.Popen(
        [samtools, "view", "-F", "2308", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert p.stdout is not None

    # Multiple ClinVar records may theoretically share a position. Keep all
    # candidate REF/ALT definitions at each position.
    candidates_at_pos = defaultdict(list)
    for row in candidates.itertuples(index=False):
        candidates_at_pos[(str(row.Chromosome), int(row.Position))].append(
            (str(row.REF), str(row.ALT))
        )

    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue

        chrom = f[2]
        targets = by_chrom.get(chrom)
        if not targets:
            continue

        try:
            mapq = int(f[4])
            start = int(f[3])
        except ValueError:
            continue

        if mapq < min_mapq:
            continue

        cigar = f[5]
        seq = f[9].upper()
        qual = f[10]

        if cigar == "*" or seq == "*" or qual == "*":
            continue

        hits = base_at_reference_positions(
            start,
            cigar,
            seq,
            qual,
            targets,
            min_bq,
        )

        for pos, base, bq in hits:
            for ref, alt in candidates_at_pos[(chrom, pos)]:
                key = (chrom, pos, ref, alt)

                if base == ref:
                    counts[key]["REF"] += 1
                elif base == alt:
                    counts[key]["ALT"] += 1
                else:
                    counts[key]["OTHER"] += 1

                bq_sum[key] += bq
                mapq_sum[key] += mapq

    rc = p.wait()
    if rc != 0:
        raise RuntimeError(
            f"samtools view failed during ClinVar base interrogation: {rc}"
        )

    rows = []

    for row in candidates.itertuples(index=False):
        key = (
            str(row.Chromosome),
            int(row.Position),
            str(row.REF),
            str(row.ALT),
        )

        ref_n = int(counts[key].get("REF", 0))
        alt_n = int(counts[key].get("ALT", 0))
        other_n = int(counts[key].get("OTHER", 0))
        total = ref_n + alt_n + other_n

        out = row._asdict()
        out.update(
            {
                "High_quality_base_observations": total,
                "REF_observations": ref_n,
                "ALT_observations": alt_n,
                "Other_observations": other_n,
                "Mean_BQ": (
                    bq_sum[key] / total if total else float("nan")
                ),
                "Mean_MAPQ": (
                    mapq_sum[key] / total if total else float("nan")
                ),
                "ALT_actually_observed": bool(alt_n > 0),
                "P_LP_ALT_observed": bool(
                    alt_n > 0 and bool(row.P_LP_annotation)
                ),
            }
        )
        rows.append(out)

    return pd.DataFrame(rows)



def _resource_first(cfg: dict, keys: list[str]) -> str | None:
    resources = cfg.get("resources", {}) or {}
    for key in keys:
        value = resources.get(key)
        if value:
            return str(value)
    return None


def _autosomal_chrom(chrom: str) -> bool:
    c = str(chrom)
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.isdigit() and 1 <= int(c) <= 22


def _filter_territory_scope(intervals, scope: str):
    if scope in {"all", "all_chromosomes", "all_nuclear_plus_mt"}:
        return list(intervals)

    if scope in {"autosomes", "autosome", "autosomal"}:
        return [x for x in intervals if _autosomal_chrom(x[0])]

    raise RuntimeError(
        "Unsupported clinvar.chromosome_scope: "
        f"{scope!r}. Use 'all' or 'autosomes'."
    )


def _sam_tag(fields, key: str):
    prefix = key + ":"

    for x in fields[11:]:
        if not x.startswith(prefix):
            continue

        parts = x.split(":", 2)
        if len(parts) != 3:
            return None

        value = parts[2]

        try:
            return int(value)
        except Exception:
            try:
                return float(value)
            except Exception:
                return value

    return None


def _aligned_query_bases(cigar: str) -> int:
    total = 0

    for ns, op in CIGAR_RE.findall(cigar):
        if op in {"M", "=", "X", "I"}:
            total += int(ns)

    return total


def _collect_plp_alt_supporting_reads(
    *,
    bam: Path,
    observed: pd.DataFrame,
    samtools: str,
    min_mapq: int,
    min_bq: int,
):
    """Recover the exact physical mates carrying candidate P/LP ALT bases."""
    if observed.empty:
        return pd.DataFrame()

    plp = observed[
        observed["P_LP_annotation"].fillna(False)
        & (pd.to_numeric(observed["ALT_observations"]) > 0)
    ].copy()

    if plp.empty:
        return pd.DataFrame()

    positions = defaultdict(list)
    variants_at_pos = defaultdict(list)

    for row in plp.itertuples(index=False):
        chrom = str(row.Chromosome)
        pos = int(row.Position)
        positions[chrom].append(pos)
        variants_at_pos[(chrom, pos)].append(
            {
                "REF": str(row.REF),
                "ALT": str(row.ALT),
                "ID": str(row.ID),
                "CLNSIG": str(row.CLNSIG),
                "CLNREVSTAT": str(row.CLNREVSTAT),
                "CLNDN": str(row.CLNDN),
                "GENEINFO": str(row.GENEINFO),
            }
        )

    for chrom in positions:
        positions[chrom] = sorted(set(positions[chrom]))

    support_rows = []
    read_ids = {}
    read_payload = {}
    next_id = 1

    p = subprocess.Popen(
        [samtools, "view", "-F", "2308", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert p.stdout is not None

    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue

        chrom = f[2]
        targets = positions.get(chrom)
        if not targets:
            continue

        try:
            flag = int(f[1])
            start = int(f[3])
            mapq = int(f[4])
        except ValueError:
            continue

        if mapq < min_mapq:
            continue

        cigar = f[5]
        seq = f[9].upper()
        qual = f[10]

        if cigar == "*" or seq == "*" or qual == "*":
            continue

        hits = base_at_reference_positions(
            start,
            cigar,
            seq,
            qual,
            targets,
            min_bq,
        )

        if not hits:
            continue

        qname = f[0]
        mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0

        read_key = (qname, mate)
        if read_key not in read_ids:
            synthetic = f"sgplp{next_id:08d}"
            next_id += 1
            read_ids[read_key] = synthetic
            read_payload[synthetic] = {
                "Sequence": seq,
                "Quality": qual,
                "QNAME": qname,
                "Mate": mate,
            }

        synthetic = read_ids[read_key]

        AS = _sam_tag(f, "AS")
        NM = _sam_tag(f, "NM")
        s1 = _sam_tag(f, "s1")
        s2 = _sam_tag(f, "s2")

        aq = _aligned_query_bases(cigar)
        identity = (
            100.0 * (1.0 - NM / aq)
            if isinstance(NM, (int, float)) and aq > 0
            else None
        )

        for pos, base, bq in hits:
            for variant in variants_at_pos[(chrom, pos)]:
                if base != variant["ALT"]:
                    continue

                support_rows.append(
                    {
                        "Synthetic_read_id": synthetic,
                        "QNAME": qname,
                        "Mate": mate,
                        "Chromosome": chrom,
                        "Position": int(pos),
                        "REF": variant["REF"],
                        "ALT": variant["ALT"],
                        "ClinVar_ID": variant["ID"],
                        "Original_alignment_start": start,
                        "Original_MAPQ": mapq,
                        "Original_CIGAR": cigar,
                        "Original_ALT_BQ": bq,
                        "Original_AS": AS,
                        "Original_NM": NM,
                        "Original_identity_percent": identity,
                        "Original_s1": s1,
                        "Original_s2": s2,
                        "Original_s2_over_s1": (
                            s2 / s1
                            if (
                                isinstance(s1, (int, float))
                                and isinstance(s2, (int, float))
                                and s1 > 0
                            )
                            else None
                        ),
                        "CLNSIG": variant["CLNSIG"],
                        "CLNREVSTAT": variant["CLNREVSTAT"],
                        "CLNDN": variant["CLNDN"],
                        "GENEINFO": variant["GENEINFO"],
                    }
                )

    rc = p.wait()
    if rc != 0:
        raise RuntimeError(
            "samtools view failed while collecting P/LP supporting reads: "
            f"{rc}"
        )

    return pd.DataFrame(support_rows), read_payload


def _write_support_fastq(read_payload: dict, path: Path):
    with path.open("w") as fh:
        for synthetic in sorted(read_payload):
            x = read_payload[synthetic]
            fh.write(f"@{synthetic}\n")
            fh.write(x["Sequence"] + "\n")
            fh.write("+\n")
            fh.write(x["Quality"] + "\n")


def _remap_supporting_reads(
    *,
    fastq: Path,
    human_mmi: Path,
    minimap2: str,
    threads: int,
    secondary_n: int,
    secondary_p: float,
):
    cmd = [
        minimap2,
        "-a",
        "-x",
        "sr",
        "-N",
        str(secondary_n),
        "-p",
        str(secondary_p),
        "--secondary=yes",
        "-t",
        str(threads),
        str(human_mmi),
        str(fastq),
    ]

    p = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )

    maps = defaultdict(list)

    for line in p.stdout.splitlines():
        if not line or line.startswith("@"):
            continue

        f = line.split("\t")
        if len(f) < 11:
            continue

        try:
            flag = int(f[1])
        except ValueError:
            continue

        if flag & 0x4:
            continue

        try:
            start = int(f[3])
            mapq = int(f[4])
        except ValueError:
            continue

        cigar = f[5]
        end = _reference_end(start, cigar)

        maps[f[0]].append(
            {
                "Flag": flag,
                "Secondary": bool(flag & 0x100),
                "Supplementary": bool(flag & 0x800),
                "Chromosome": f[2],
                "Start": start,
                "End": end,
                "MAPQ": mapq,
                "CIGAR": cigar,
                "Sequence": f[9],
                "Quality": f[10],
                "AS": _sam_tag(f, "AS"),
                "NM": _sam_tag(f, "NM"),
                "s1": _sam_tag(f, "s1"),
                "s2": _sam_tag(f, "s2"),
            }
        )

    return maps


def _score_value(x):
    value = x.get("AS")
    if isinstance(value, (int, float)):
        return float(value)
    return float("-inf")


def _alignment_alt_at_target(aln: dict, target_pos: int, alt: str):
    if aln.get("Sequence") in {None, "*"}:
        return None
    if aln.get("Quality") in {None, "*"}:
        return None

    hits = base_at_reference_positions(
        int(aln["Start"]),
        str(aln["CIGAR"]),
        str(aln["Sequence"]).upper(),
        str(aln["Quality"]),
        [int(target_pos)],
        0,
    )

    if not hits:
        return False

    return any(base == alt for _, base, _ in hits)


def _classify_plp_support_rows(
    support: pd.DataFrame,
    maps: dict,
    *,
    min_remap_mapq: int,
    distinct_locus_window_bp: int,
):
    if support.empty:
        return support.copy()

    rows = []

    for x in support.to_dict(orient="records"):
        candidates = list(maps.get(x["Synthetic_read_id"], []))
        candidates.sort(
            key=lambda z: (_score_value(z), z.get("MAPQ", 0)),
            reverse=True,
        )

        target = [
            z
            for z in candidates
            if (
                z["Chromosome"] == x["Chromosome"]
                and int(z["Start"]) <= int(x["Position"]) <= int(z["End"])
            )
        ]
        target.sort(
            key=lambda z: (_score_value(z), z.get("MAPQ", 0)),
            reverse=True,
        )

        best_target = target[0] if target else None
        best_overall = candidates[0] if candidates else None

        competitor = None
        if best_target is not None:
            for z in candidates:
                same_local_locus = (
                    z["Chromosome"] == best_target["Chromosome"]
                    and abs(
                        int(z["Start"]) - int(best_target["Start"])
                    )
                    <= distinct_locus_window_bp
                )
                if not same_local_locus:
                    competitor = z
                    break

        target_as = (
            _score_value(best_target)
            if best_target is not None
            else float("-inf")
        )
        overall_as = (
            _score_value(best_overall)
            if best_overall is not None
            else float("-inf")
        )
        competitor_as = (
            _score_value(competitor)
            if competitor is not None
            else None
        )

        target_is_top_score = (
            best_target is not None and target_as == overall_as
        )
        no_equal_or_better_distinct_competitor = (
            best_target is not None
            and (
                competitor is None
                or competitor_as < target_as
            )
        )
        alt_preserved = (
            _alignment_alt_at_target(
                best_target,
                int(x["Position"]),
                str(x["ALT"]),
            )
            if best_target is not None
            else False
        )

        pass_validation = bool(
            best_target is not None
            and alt_preserved is True
            and int(best_target["MAPQ"]) >= min_remap_mapq
            and target_is_top_score
            and no_equal_or_better_distinct_competitor
        )

        ratio = None
        if (
            best_target is not None
            and competitor is not None
            and isinstance(target_as, (int, float))
            and target_as > 0
            and isinstance(competitor_as, (int, float))
        ):
            ratio = competitor_as / target_as

        out = dict(x)
        out.update(
            {
                "Remap_alignments": len(candidates),
                "Target_alignment_found": bool(best_target is not None),
                "Target_ALT_preserved_on_remap": alt_preserved,
                "Target_is_top_AS": bool(target_is_top_score),
                "Target_best_start": (
                    best_target["Start"] if best_target else None
                ),
                "Target_best_AS": (
                    best_target["AS"] if best_target else None
                ),
                "Target_best_MAPQ": (
                    best_target["MAPQ"] if best_target else None
                ),
                "Distinct_competitor_chromosome": (
                    competitor["Chromosome"] if competitor else ""
                ),
                "Distinct_competitor_start": (
                    competitor["Start"] if competitor else None
                ),
                "Distinct_competitor_AS": (
                    competitor["AS"] if competitor else None
                ),
                "Distinct_competitor_MAPQ": (
                    competitor["MAPQ"] if competitor else None
                ),
                "Competitor_over_target_AS": ratio,
                "Unambiguous_locus_support": pass_validation,
            }
        )
        rows.append(out)

    return pd.DataFrame(rows)



def _extract_supporting_pairs_from_bam(
    *,
    bam: Path,
    support: pd.DataFrame,
    samtools: str,
    out_r1: Path,
    out_r2: Path,
):
    """Extract full physical pairs for candidate P/LP supporting QNAMEs."""
    if support.empty:
        return 0

    qnames = sorted(set(support["QNAME"].astype(str)))
    if not qnames:
        return 0

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".qnames",
        dir=Path(out_r1).parent,
        delete=False,
    ) as tf:
        qfile = Path(tf.name)
        for qn in qnames:
            tf.write(qn + "\n")

    try:
        cmd = (
            f"{samtools} view -b -N {qfile} {bam} "
            f"| {samtools} collate -u -O - "
            f"| {samtools} fastq -n "
            f"-1 {out_r1} -2 {out_r2} "
            f"-0 /dev/null -s /dev/null -"
        )

        p = subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            text=True,
            capture_output=True,
        )

        if p.returncode != 0:
            raise RuntimeError(
                "Paired FASTQ extraction failed:\n"
                + p.stderr[-4000:]
            )

    finally:
        qfile.unlink(missing_ok=True)

    return len(qnames)


def _remap_pairs(
    *,
    r1: Path,
    r2: Path,
    human_mmi: Path,
    minimap2: str,
    threads: int,
    secondary_n: int,
    secondary_p: float,
):
    """Remap supporting physical pairs genome-wide with secondary hits."""
    cmd = [
        minimap2,
        "-a",
        "-x",
        "sr",
        "-N",
        str(secondary_n),
        "-p",
        str(secondary_p),
        "--secondary=yes",
        "-t",
        str(threads),
        str(human_mmi),
        str(r1),
        str(r2),
    ]

    p = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )

    maps = defaultdict(list)

    for line in p.stdout.splitlines():
        if not line or line.startswith("@"):
            continue

        f = line.split("\t")
        if len(f) < 11:
            continue

        try:
            flag = int(f[1])
            start = int(f[3])
            mapq = int(f[4])
        except ValueError:
            continue

        if flag & 0x4:
            continue

        cigar = f[5]
        end = _reference_end(start, cigar)

        mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0

        maps[f[0]].append(
            {
                "Flag": flag,
                "Mate": mate,
                "Secondary": bool(flag & 0x100),
                "Supplementary": bool(flag & 0x800),
                "Proper_pair": bool(flag & 0x2),
                "Chromosome": f[2],
                "Start": start,
                "End": end,
                "MAPQ": mapq,
                "CIGAR": cigar,
                "Sequence": f[9],
                "Quality": f[10],
                "AS": _sam_tag(f, "AS"),
                "NM": _sam_tag(f, "NM"),
                "s1": _sam_tag(f, "s1"),
                "s2": _sam_tag(f, "s2"),
            }
        )

    return maps


def _classify_pair_support_rows(
    support: pd.DataFrame,
    pair_maps: dict,
    *,
    min_remap_mapq: int,
    distinct_locus_window_bp: int,
):
    """Evaluate candidate ALT support using the full physical read pair.

    Canonical paired-end support requires the ALT-supporting mate to map back
    over the target locus with the ALT preserved and with sufficient MAPQ.
    The pair must be primary/proper at that placement. Distinct equal-or-better
    target-mate competitors are rejected when present.

    Independent-mate validation remains a separate, stronger diagnostic.
    """
    if support.empty:
        return support.copy()

    rows = []

    for x in support.to_dict(orient="records"):
        qname = str(x["QNAME"])
        target_mate = int(x["Mate"])
        candidates = list(pair_maps.get(qname, []))

        # Focus on the same physical mate that originally carried the ALT.
        mate_maps = [
            z for z in candidates
            if int(z.get("Mate", 0)) == target_mate
        ]

        mate_maps.sort(
            key=lambda z: (_score_value(z), z.get("MAPQ", 0)),
            reverse=True,
        )

        target = [
            z
            for z in mate_maps
            if (
                z["Chromosome"] == x["Chromosome"]
                and int(z["Start"]) <= int(x["Position"]) <= int(z["End"])
            )
        ]
        target.sort(
            key=lambda z: (_score_value(z), z.get("MAPQ", 0)),
            reverse=True,
        )

        best_target = target[0] if target else None
        best_overall = mate_maps[0] if mate_maps else None

        competitor = None

        if best_target is not None:
            for z in mate_maps:
                same_local_locus = (
                    z["Chromosome"] == best_target["Chromosome"]
                    and abs(
                        int(z["Start"]) - int(best_target["Start"])
                    )
                    <= distinct_locus_window_bp
                )

                if not same_local_locus:
                    competitor = z
                    break

        target_as = (
            _score_value(best_target)
            if best_target is not None
            else float("-inf")
        )

        overall_as = (
            _score_value(best_overall)
            if best_overall is not None
            else float("-inf")
        )

        competitor_as = (
            _score_value(competitor)
            if competitor is not None
            else None
        )

        target_is_top_score = (
            best_target is not None
            and target_as == overall_as
        )

        no_equal_or_better_distinct_competitor = (
            best_target is not None
            and (
                competitor is None
                or competitor_as < target_as
            )
        )

        alt_preserved = (
            _alignment_alt_at_target(
                best_target,
                int(x["Position"]),
                str(x["ALT"]),
            )
            if best_target is not None
            else False
        )

        canonical_pair_support = bool(
            best_target is not None
            and bool(best_target.get("Proper_pair", False))
            and alt_preserved is True
            and int(best_target["MAPQ"]) >= min_remap_mapq
            and target_is_top_score
            and no_equal_or_better_distinct_competitor
        )

        ratio = None

        if (
            best_target is not None
            and competitor is not None
            and isinstance(target_as, (int, float))
            and target_as > 0
            and isinstance(competitor_as, (int, float))
        ):
            ratio = competitor_as / target_as

        out = dict(x)
        out.update(
            {
                "Pair_remap_alignments": len(candidates),
                "Pair_target_alignment_found": bool(
                    best_target is not None
                ),
                "Pair_target_ALT_preserved": alt_preserved,
                "Pair_target_is_top_AS": bool(target_is_top_score),
                "Pair_target_proper_pair": (
                    bool(best_target.get("Proper_pair", False))
                    if best_target is not None
                    else False
                ),
                "Pair_target_best_start": (
                    best_target["Start"] if best_target else None
                ),
                "Pair_target_best_AS": (
                    best_target["AS"] if best_target else None
                ),
                "Pair_target_best_MAPQ": (
                    best_target["MAPQ"] if best_target else None
                ),
                "Pair_distinct_competitor_chromosome": (
                    competitor["Chromosome"] if competitor else ""
                ),
                "Pair_distinct_competitor_start": (
                    competitor["Start"] if competitor else None
                ),
                "Pair_distinct_competitor_AS": (
                    competitor["AS"] if competitor else None
                ),
                "Pair_competitor_over_target_AS": ratio,
                "Pair_locus_supported": canonical_pair_support,
            }
        )

        rows.append(out)

    return pd.DataFrame(rows)


def _merge_independent_and_pair_validation(
    independent: pd.DataFrame,
    pair_validation: pd.DataFrame,
):
    """Combine stronger independent-mate diagnostics with canonical pair QC."""
    if pair_validation.empty:
        return independent.copy()

    keys = [
        "Chromosome",
        "Position",
        "REF",
        "ALT",
        "QNAME",
        "Mate",
    ]

    if independent.empty:
        out = pair_validation.copy()
        out["Independent_mate_locus_supported"] = False
        return out

    independent = independent.copy().rename(
        columns={
            "Unambiguous_locus_support":
                "Independent_mate_locus_supported",
        }
    )

    keep_pair = keys + [
        c for c in pair_validation.columns
        if c.startswith("Pair_")
    ]

    return independent.merge(
        pair_validation[keep_pair],
        on=keys,
        how="outer",
    )


def _attach_site_validation(
    observed: pd.DataFrame,
    audit: pd.DataFrame,
    *,
    min_passing_mates: int,
):
    out = observed.copy()

    out["Candidate_P_LP_ALT_observed"] = (
        out["P_LP_annotation"].fillna(False)
        & (pd.to_numeric(out["ALT_observations"]) > 0)
    )

    out["P_LP_ALT_supporting_mates"] = 0
    out["P_LP_pair_locus_supported_mates"] = 0
    out["P_LP_independent_mate_supported_mates"] = 0
    out["Validated_P_LP_ALT"] = False

    if audit.empty:
        return out

    keys = ["Chromosome", "Position", "REF", "ALT"]

    working = audit.copy()

    if "Pair_locus_supported" not in working.columns:
        working["Pair_locus_supported"] = False

    if "Independent_mate_locus_supported" not in working.columns:
        if "Unambiguous_locus_support" in working.columns:
            working["Independent_mate_locus_supported"] = working[
                "Unambiguous_locus_support"
            ].fillna(False)
        else:
            working["Independent_mate_locus_supported"] = False

    grouped = (
        working.groupby(keys, dropna=False)
        .agg(
            P_LP_ALT_supporting_mates=(
                "QNAME",
                "nunique",
            ),
            P_LP_pair_locus_supported_mates=(
                "Pair_locus_supported",
                "sum",
            ),
            P_LP_independent_mate_supported_mates=(
                "Independent_mate_locus_supported",
                "sum",
            ),
        )
        .reset_index()
    )

    out = out.drop(
        columns=[
            "P_LP_ALT_supporting_mates",
            "P_LP_pair_locus_supported_mates",
            "P_LP_independent_mate_supported_mates",
        ],
        errors="ignore",
    ).merge(
        grouped,
        on=keys,
        how="left",
    )

    for col in (
        "P_LP_ALT_supporting_mates",
        "P_LP_pair_locus_supported_mates",
        "P_LP_independent_mate_supported_mates",
    ):
        out[col] = pd.to_numeric(
            out[col],
            errors="coerce",
        ).fillna(0).astype(int)

    out["Validated_P_LP_ALT"] = (
        out["Candidate_P_LP_ALT_observed"]
        & (
            out["P_LP_pair_locus_supported_mates"]
            >= int(min_passing_mates)
        )
    )

    return out


def _empty_summary(sample: str):
    return {
        "Sample": sample,
        "Strict_fragment_territory_intervals": 0,
        "Candidate_ClinVar_biallelic_SNVs": 0,
        "ClinVar_sites_with_HQ_base_observation": 0,
        "ClinVar_ALT_sites_observed": 0,
        "ClinVar_ALT_base_observations": 0,
        "ClinVar_candidate_P_LP_ALT_sites": 0,
        "ClinVar_candidate_P_LP_ALT_base_observations": 0,
        "ClinVar_validated_P_LP_ALT_sites": 0,
        "ClinVar_validated_P_LP_supporting_mates": 0,
        "ClinVar_exposure_detected": False,
        "P_LP_candidate_exposure_detected": False,
        "P_LP_validated_exposure_detected": False,
    }



def _count_exposure(observed: pd.DataFrame) -> dict:
    if observed.empty:
        return {
            "sites_with_hq": 0,
            "alt_sites": 0,
            "alt_bases": 0,
            "candidate_plp_sites": 0,
            "candidate_plp_alt_bases": 0,
        }

    hq = pd.to_numeric(
        observed.get(
            "High_quality_base_observations",
            pd.Series(dtype=int),
        ),
        errors="coerce",
    ).fillna(0)

    alt = pd.to_numeric(
        observed.get(
            "ALT_observations",
            pd.Series(dtype=int),
        ),
        errors="coerce",
    ).fillna(0)

    candidate_plp = (
        observed.get(
            "P_LP_annotation",
            pd.Series(False, index=observed.index),
        ).fillna(False)
        & (alt > 0)
    )

    return {
        "sites_with_hq": int((hq > 0).sum()),
        "alt_sites": int((alt > 0).sum()),
        "alt_bases": int(alt.sum()),
        "candidate_plp_sites": int(candidate_plp.sum()),
        "candidate_plp_alt_bases": int(
            alt.loc[candidate_plp].sum()
            if candidate_plp.any()
            else 0
        ),
    }


def _resolve_resource_path(value: str | None) -> Path | None:
    if not value:
        return None

    p = Path(value)
    if p.is_absolute():
        return p

    # Current repository root is the most natural resolution for bundled
    # resources; package-root fallback keeps invocation outside the repo sane.
    candidates = [
        Path.cwd() / p,
        Path(__file__).resolve().parents[2] / p,
    ]

    for x in candidates:
        if x.exists():
            return x.resolve()

    return candidates[0].resolve()


def run_clinvar(
    *,
    output: str | Path,
    sample: str,
    config_path: str | None = None,
) -> dict:
    output = Path(output).resolve()
    cfg, resolved_cfg = load_effective_config(config_path)
    settings = _settings(cfg)

    if not settings["enabled"]:
        return {"enabled": False}

    vcf_value = _cfg_get(cfg, "resources", "clinvar_vcf", default=None)
    hotspot_value = _resource_first(
        cfg,
        [
            "hotspot_exclusion_bed",
            "hotspot_exclusion_mask",
        ],
    )

    if not vcf_value:
        raise RuntimeError("resources.clinvar_vcf is not configured")

    vcf = Path(vcf_value)
    if not vcf.is_file():
        raise FileNotFoundError(f"ClinVar VCF missing: {vcf}")

    evidence = output / "evidence"
    strict_bam = evidence / "strict_human.bam"

    if not strict_bam.is_file():
        raise FileNotFoundError(
            f"Retained strict-human BAM not found: {strict_bam}"
        )

    samtools = _tool(cfg, "samtools", "samtools")
    bcftools = _tool(cfg, "bcftools", "bcftools")

    territory_raw = _strict_fragment_territory(
        bam=strict_bam,
        samtools=samtools,
        min_mapq=settings["min_mapq"],
    )
    territory_raw = _filter_territory_scope(
        territory_raw,
        settings["chromosome_scope"],
    )

    hotspot = []
    hotspot_path = None

    if settings["apply_hotspot_mask"]:
        hotspot_path = _resolve_resource_path(hotspot_value)

        if hotspot_path is None:
            raise RuntimeError(
                "ClinVar hotspot exclusion is enabled but neither "
                "resources.hotspot_exclusion_bed nor "
                "resources.hotspot_exclusion_mask is configured."
            )

        if not hotspot_path.is_file():
            raise FileNotFoundError(
                f"Hotspot exclusion mask missing: {hotspot_path}"
            )

        hotspot = _read_bed(hotspot_path)

    territory_filtered = (
        _subtract_bed(territory_raw, hotspot)
        if settings["apply_hotspot_mask"]
        else list(territory_raw)
    )

    runtime_root_value = _cfg_get(
        cfg, "runtime", "tempdir", default=None
    )
    runtime_root = (
        Path(runtime_root_value)
        if runtime_root_value
        else Path(tempfile.gettempdir()) / "srafeguard"
    )
    runtime_root.mkdir(parents=True, exist_ok=True)

    raw_site_out = output / "clinvar_site_observations_raw.tsv.gz"
    site_out = output / "clinvar_site_observations.tsv.gz"
    audit_out = output / "clinvar_plp_read_validation.tsv.gz"
    summary_out = output / "clinvar_exposure.tsv"
    metadata_out = output / "clinvar_analysis_metadata.json"

    raw_candidates = pd.DataFrame()
    raw_observed = pd.DataFrame()
    filtered_candidates = pd.DataFrame()
    filtered_observed = pd.DataFrame()
    audit = pd.DataFrame()

    with tempfile.TemporaryDirectory(
        prefix=f"srafeguard_clinvar_{sample}_",
        dir=runtime_root,
    ) as td:
        td = Path(td)
        from .owned_temporary import register
        register(output,td,"clinical_workspace")

        # --------------------------------------------------------------
        # RAW exposure: preserve what is actually recoverable before
        # recurrent-hotspot QC.
        # --------------------------------------------------------------
        if territory_raw:
            raw_bed = td / "strict_fragment_territory_raw.bed"
            _write_bed(territory_raw, raw_bed)

            raw_candidates = _query_clinvar(
                vcf=vcf,
                territory_bed=raw_bed,
                bcftools=bcftools,
                biallelic_snv_only=settings["biallelic_snv_only"],
            )

            raw_observed = _observe_candidate_sites(
                bam=strict_bam,
                candidates=raw_candidates,
                samtools=samtools,
                min_mapq=settings["min_mapq"],
                min_bq=settings["min_bq"],
            )

        # --------------------------------------------------------------
        # CANONICAL high-confidence exposure after recurrent-hotspot QC.
        # Generic locus-specific P/LP validation runs only here.
        # --------------------------------------------------------------
        if territory_filtered:
            filtered_bed = td / "strict_fragment_territory_post_hotspot.bed"
            _write_bed(territory_filtered, filtered_bed)

            filtered_candidates = _query_clinvar(
                vcf=vcf,
                territory_bed=filtered_bed,
                bcftools=bcftools,
                biallelic_snv_only=settings["biallelic_snv_only"],
            )

            filtered_observed = _observe_candidate_sites(
                bam=strict_bam,
                candidates=filtered_candidates,
                samtools=samtools,
                min_mapq=settings["min_mapq"],
                min_bq=settings["min_bq"],
            )

            candidate_mask = (
                filtered_observed["P_LP_annotation"].fillna(False)
                & (
                    pd.to_numeric(
                        filtered_observed["ALT_observations"],
                        errors="coerce",
                    ).fillna(0)
                    > 0
                )
            ) if not filtered_observed.empty else pd.Series(dtype=bool)

            if (
                settings["validate_plp"]
                and not filtered_observed.empty
                and candidate_mask.any()
            ):
                human_mmi_value = _resource_first(
                    cfg,
                    [
                        "human_mmi",
                        "human_minimap2_index",
                        "human_minimap_index",
                        "minimap2_human_index",
                        "human_reference_mmi",
                    ],
                )

                if not human_mmi_value:
                    raise RuntimeError(
                        "P/LP validation requires the human minimap2 index."
                    )

                human_mmi = Path(human_mmi_value)
                if not human_mmi.is_file():
                    raise FileNotFoundError(
                        f"Human minimap2 index missing: {human_mmi}"
                    )

                minimap2 = _tool(cfg, "minimap2", "minimap2")

                support, read_payload = _collect_plp_alt_supporting_reads(
                    bam=strict_bam,
                    observed=filtered_observed,
                    samtools=samtools,
                    min_mapq=settings["min_mapq"],
                    min_bq=settings["min_bq"],
                )

                if not support.empty:
                    support_fq = td / "plp_alt_supporting_mates.fastq"
                    _write_support_fastq(read_payload, support_fq)

                    remaps = _remap_supporting_reads(
                        fastq=support_fq,
                        human_mmi=human_mmi,
                        minimap2=minimap2,
                        threads=settings["validation_threads"],
                        secondary_n=settings["validation_secondary_n"],
                        secondary_p=settings["validation_secondary_p"],
                    )

                    independent_audit = _classify_plp_support_rows(
                        support,
                        remaps,
                        min_remap_mapq=settings[
                            "validation_min_remap_mapq"
                        ],
                        distinct_locus_window_bp=settings[
                            "validation_distinct_locus_window_bp"
                        ],
                    )

                    pair_r1 = td / "plp_support_pairs_R1.fastq"
                    pair_r2 = td / "plp_support_pairs_R2.fastq"

                    _extract_supporting_pairs_from_bam(
                        bam=strict_bam,
                        support=support,
                        samtools=samtools,
                        out_r1=pair_r1,
                        out_r2=pair_r2,
                    )

                    pair_maps = _remap_pairs(
                        r1=pair_r1,
                        r2=pair_r2,
                        human_mmi=human_mmi,
                        minimap2=minimap2,
                        threads=settings["validation_threads"],
                        secondary_n=settings["validation_secondary_n"],
                        secondary_p=settings["validation_secondary_p"],
                    )

                    pair_audit = _classify_pair_support_rows(
                        support,
                        pair_maps,
                        min_remap_mapq=settings[
                            "validation_min_remap_mapq"
                        ],
                        distinct_locus_window_bp=settings[
                            "validation_distinct_locus_window_bp"
                        ],
                    )

                    audit = _merge_independent_and_pair_validation(
                        independent_audit,
                        pair_audit,
                    )

    # Attach locus-validation status only to the canonical post-hotspot table.
    if filtered_observed.empty:
        filtered_observed = pd.DataFrame()
    else:
        filtered_observed = _attach_site_validation(
            filtered_observed,
            audit,
            min_passing_mates=settings[
                "validation_min_passing_mates"
            ],
        )

    raw_observed.to_csv(
        raw_site_out,
        sep="\t",
        index=False,
        compression="gzip",
    )

    filtered_observed.to_csv(
        site_out,
        sep="\t",
        index=False,
        compression="gzip",
    )

    audit.to_csv(
        audit_out,
        sep="\t",
        index=False,
        compression="gzip",
    )

    raw = _count_exposure(raw_observed)
    post = _count_exposure(filtered_observed)

    hotspot_alt_sites = max(0, raw["alt_sites"] - post["alt_sites"])
    hotspot_alt_bases = max(0, raw["alt_bases"] - post["alt_bases"])
    hotspot_plp_sites = max(
        0,
        raw["candidate_plp_sites"] - post["candidate_plp_sites"],
    )
    hotspot_plp_bases = max(
        0,
        raw["candidate_plp_alt_bases"]
        - post["candidate_plp_alt_bases"],
    )

    validated_plp = (
        filtered_observed.get(
            "Validated_P_LP_ALT",
            pd.Series(dtype=bool),
        ).fillna(False)
        if not filtered_observed.empty
        else pd.Series(dtype=bool)
    )

    validated_support = (
        int(
            filtered_observed.loc[
                validated_plp,
                "P_LP_pair_locus_supported_mates",
            ].sum()
        )
        if (
            not filtered_observed.empty
            and validated_plp.any()
            and "P_LP_pair_locus_supported_mates"
            in filtered_observed.columns
        )
        else 0
    )

    summary = {
        "Sample": sample,

        # Raw/privacy-recoverability layer
        "Raw_Strict_fragment_territory_intervals": int(len(territory_raw)),
        "Raw_Candidate_ClinVar_biallelic_SNVs": int(len(raw_candidates)),
        "Raw_ClinVar_sites_with_HQ_base_observation": raw["sites_with_hq"],
        "Raw_ClinVar_ALT_sites_observed": raw["alt_sites"],
        "Raw_ClinVar_ALT_base_observations": raw["alt_bases"],
        "Raw_ClinVar_candidate_P_LP_ALT_sites": raw["candidate_plp_sites"],
        "Raw_ClinVar_candidate_P_LP_ALT_base_observations": (
            raw["candidate_plp_alt_bases"]
        ),

        # Explicit hotspot layer
        "Hotspot_excluded_territory_intervals": int(
            max(0, len(territory_raw) - len(territory_filtered))
        ),
        "Hotspot_ClinVar_ALT_sites_observed": hotspot_alt_sites,
        "Hotspot_ClinVar_ALT_base_observations": hotspot_alt_bases,
        "Hotspot_candidate_P_LP_ALT_sites": hotspot_plp_sites,
        "Hotspot_candidate_P_LP_ALT_base_observations": hotspot_plp_bases,

        # Canonical post-QC layer
        "Strict_fragment_territory_intervals": int(len(territory_filtered)),
        "Candidate_ClinVar_biallelic_SNVs": int(len(filtered_candidates)),
        "ClinVar_sites_with_HQ_base_observation": post["sites_with_hq"],
        "ClinVar_ALT_sites_observed": post["alt_sites"],
        "ClinVar_ALT_base_observations": post["alt_bases"],
        "ClinVar_candidate_P_LP_ALT_sites": post["candidate_plp_sites"],
        "ClinVar_candidate_P_LP_ALT_base_observations": (
            post["candidate_plp_alt_bases"]
        ),
        "ClinVar_validated_P_LP_ALT_sites": int(validated_plp.sum()),
        "ClinVar_validated_P_LP_supporting_mates": validated_support,

        "Raw_ClinVar_exposure_detected": bool(raw["alt_sites"] > 0),
        "ClinVar_exposure_detected": bool(post["alt_sites"] > 0),
        "P_LP_candidate_exposure_detected": bool(
            post["candidate_plp_sites"] > 0
        ),
        "P_LP_validated_exposure_detected": bool(validated_plp.any()),
    }

    pd.DataFrame([summary]).to_csv(
        summary_out,
        sep="\t",
        index=False,
    )

    metadata = {
        "sample": sample,
        "config": resolved_cfg,
        "chromosome_scope": settings["chromosome_scope"],
        "clinvar_vcf": str(vcf),
        "apply_hotspot_exclusion_mask": settings["apply_hotspot_mask"],
        "hotspot_exclusion_bed": (
            str(hotspot_path) if hotspot_path else None
        ),
        "territory_intervals_raw": len(territory_raw),
        "territory_intervals_post_hotspot": len(territory_filtered),
        "min_mapq": settings["min_mapq"],
        "min_bq": settings["min_bq"],
        "biallelic_snv_only": settings["biallelic_snv_only"],
        "raw_candidate_clinvar_snvs": int(len(raw_candidates)),
        "post_hotspot_candidate_clinvar_snvs": int(
            len(filtered_candidates)
        ),
        "raw_clinvar_alt_sites": raw["alt_sites"],
        "post_hotspot_clinvar_alt_sites": post["alt_sites"],
        "raw_candidate_plp_alt_sites": raw["candidate_plp_sites"],
        "hotspot_candidate_plp_alt_sites": hotspot_plp_sites,
        "post_hotspot_candidate_plp_alt_sites": post[
            "candidate_plp_sites"
        ],
        "validated_plp_alt_sites": int(validated_plp.sum()),
        "validation": {
            "enabled": settings["validate_plp"],
            "secondary_n": settings["validation_secondary_n"],
            "secondary_p": settings["validation_secondary_p"],
            "min_remap_mapq": settings[
                "validation_min_remap_mapq"
            ],
            "distinct_locus_window_bp": settings[
                "validation_distinct_locus_window_bp"
            ],
            "min_passing_mates": settings[
                "validation_min_passing_mates"
            ],
            "rule": (
                "For paired-end input, post-hotspot candidate P/LP sites are "
                "canonically validated using the full physical read pair. "
                "The ALT-supporting mate must remap over the target locus "
                "with the ALT preserved, sufficient MAPQ, proper pairing, "
                "top alignment score, and no distinct equal-or-better "
                "competitor. Independent-mate remapping is retained as a "
                "stricter diagnostic, not a required criterion."
            ),
        },
    }

    metadata_out.write_text(json.dumps(metadata, indent=2))

    _print_summary(summary)
    return summary


def _print_summary(summary: dict):
    print()
    print("=" * 80)
    print("SRAFEGUARD — CLINVAR EXPOSURE")
    print("=" * 80)
    print(
        f"Raw ClinVar ALT sites:          "
        f"{summary['Raw_ClinVar_ALT_sites_observed']:,}"
    )
    print(
        f"Raw candidate P/LP ALT sites:   "
        f"{summary['Raw_ClinVar_candidate_P_LP_ALT_sites']:,}"
    )
    print(
        f"ALT sites in excluded hotspots: "
        f"{summary['Hotspot_ClinVar_ALT_sites_observed']:,}"
    )
    print(
        f"P/LP candidates in hotspots:    "
        f"{summary['Hotspot_candidate_P_LP_ALT_sites']:,}"
    )
    print(
        f"Post-QC ClinVar ALT sites:       "
        f"{summary['ClinVar_ALT_sites_observed']:,}"
    )
    print(
        f"Post-QC candidate P/LP sites:    "
        f"{summary['ClinVar_candidate_P_LP_ALT_sites']:,}"
    )
    print(
        f"Computational P/LP placement sites:        "
        f"{summary['ClinVar_validated_P_LP_ALT_sites']:,}"
    )
    print(
        f"Computational support mates:         "
        f"{summary['ClinVar_validated_P_LP_supporting_mates']:,}"
    )
    print(
        "Raw clinical-variant exposure:   "
        + (
            "DETECTED"
            if summary["Raw_ClinVar_exposure_detected"]
            else "NOT DETECTED"
        )
    )
    print(
        "Post-QC clinical exposure:       "
        + (
            "DETECTED"
            if summary["ClinVar_exposure_detected"]
            else "NOT DETECTED"
        )
    )
    print(
        "Computational P/LP placement:         "
        + (
            "DETECTED"
            if summary["P_LP_validated_exposure_detected"]
            else "NOT DETECTED"
        )
    )
