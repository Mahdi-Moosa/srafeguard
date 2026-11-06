from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .runtime_config import load_effective_config

DEFAULT_POPS = ["AFR", "EUR", "SAS", "EAS", "AMR"]
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
        raise FileNotFoundError(f"Required tool '{fallback}' not found in PATH.")
    return found


def _settings(cfg: dict) -> dict:
    p = cfg.get("population", {}) or {}
    pops = [str(x).upper() for x in (p.get("populations") or DEFAULT_POPS)]
    return {
        "enabled": bool(p.get("enabled", True)),
        "extraction_mode": str(p.get("extraction_mode", "original_alignment")),
        "populations": pops,
        "min_pairwise_af_delta": float(p.get("min_pairwise_af_delta", 0.30)),
        "af_floor": float(p.get("af_floor", 0.005)),
        "min_mapq": int(p.get("min_mapq", 30)),
        "min_bq": int(p.get("min_bq", 30)),
        "block_size_bp": int(p.get("block_size_bp", 1_000_000)),
        "bootstraps": int(p.get("bootstraps", 1000)),
        "seed": int(p.get("seed", 20260902)),
        "informative_min_sites": int(p.get("informative_min_sites", 20)),
        "informative_min_bootstrap": float(p.get("informative_min_bootstrap", 0.90)),
    }


def _autosomal(chrom: str) -> bool:
    c = str(chrom)
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.isdigit() and 1 <= int(c) <= 22


def _marker_cache(marker_path: Path, cache_root: Path, delta: float, pops: list[str]) -> Path:
    stat = marker_path.stat()
    token = "|".join([
        str(marker_path.resolve()), str(stat.st_size), str(stat.st_mtime_ns),
        f"{delta:.8f}", ",".join(pops)
    ])
    key = hashlib.sha256(token.encode()).hexdigest()[:16]
    out = cache_root / f"population_markers_{key}.tsv.gz"
    if out.is_file():
        return out

    cache_root.mkdir(parents=True, exist_ok=True)
    markers = pd.read_csv(marker_path, sep="\t", compression="infer")
    required = {"Chromosome", "Position", "REF", "ALT", "Max_pairwise_AF_delta"}
    required |= {f"AF_{p}" for p in pops}
    missing = sorted(required - set(markers.columns))
    if missing:
        raise RuntimeError("Population marker resource is missing columns: " + ", ".join(missing))

    markers = markers[
        markers["Chromosome"].map(_autosomal)
        & (pd.to_numeric(markers["Max_pairwise_AF_delta"], errors="coerce") >= delta)
    ].copy()
    markers["Position"] = pd.to_numeric(markers["Position"], errors="coerce").astype("Int64")
    markers["REF"] = markers["REF"].astype(str).str.upper()
    markers["ALT"] = markers["ALT"].astype(str).str.upper()
    markers = markers[
        markers["Position"].notna()
        & markers["REF"].str.fullmatch("[ACGT]")
        & markers["ALT"].str.fullmatch("[ACGT]")
    ].copy()
    markers["Position"] = markers["Position"].astype(int)
    keep = ["Chromosome", "Position", "REF", "ALT", "Max_pairwise_AF_delta"] + [f"AF_{p}" for p in pops]
    markers = markers[keep].sort_values(["Chromosome", "Position"]).drop_duplicates(["Chromosome", "Position"])
    tmp = out.with_name(out.name + ".tmp")
    markers.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, out)
    return out


def base_at_reference_positions(aln_start: int, cigar: str, seq: str, qual: str,
                                target_positions: list[int], min_bq: int):
    """Frozen Phase3A-style manual CIGAR projection."""
    out = []
    rpos, qpos = aln_start, 0
    for ns, op in CIGAR_RE.findall(cigar):
        n = int(ns)
        if op in {"M", "=", "X"}:
            start, end = rpos, rpos + n - 1
            lo = bisect_left(target_positions, start)
            hi = bisect_right(target_positions, end)
            for pos in target_positions[lo:hi]:
                qi = qpos + (pos - rpos)
                if not (0 <= qi < len(seq)) or qual == "*" or qi >= len(qual):
                    continue
                base = seq[qi].upper()
                if base not in {"A", "C", "G", "T"}:
                    continue
                bq = ord(qual[qi]) - 33
                if bq >= min_bq:
                    out.append((pos, base, bq))
            rpos += n
            qpos += n
        elif op in {"I", "S"}:
            qpos += n
        elif op in {"D", "N"}:
            rpos += n
        elif op in {"H", "P"}:
            pass
    return out


def _empty_obs(pops: list[str]) -> pd.DataFrame:
    cols = [
        "Chromosome", "Position", "REF", "ALT", "Observed_state", "Observed_allele",
        "REF_observations", "ALT_observations", "Total_usable_base_observations",
        "Mean_BQ", "Mean_MAPQ", "Block_1Mb", "Max_pairwise_AF_delta"
    ] + [f"AF_{p}" for p in pops]
    return pd.DataFrame(columns=cols)


def _extract_from_original_bam(*, bam: Path, markers: pd.DataFrame, samtools: str,
                               min_mapq: int, min_bq: int, block_size_bp: int,
                               pops: list[str]) -> pd.DataFrame:
    marker_positions: dict[str, list[int]] = {}
    lookup: dict[tuple[str, int], dict] = {}
    for row in markers.itertuples(index=False):
        chrom, pos = str(row.Chromosome), int(row.Position)
        marker_positions.setdefault(chrom, []).append(pos)
        lookup[(chrom, pos)] = {
            "REF": str(row.REF).upper(), "ALT": str(row.ALT).upper(),
            "Max_pairwise_AF_delta": float(row.Max_pairwise_AF_delta),
            **{f"AF_{p}": float(getattr(row, f"AF_{p}")) for p in pops},
        }
    for chrom in marker_positions:
        marker_positions[chrom] = sorted(set(marker_positions[chrom]))

    obs = defaultdict(lambda: defaultdict(int))
    bq_sum, mapq_sum = defaultdict(int), defaultdict(int)

    p = subprocess.Popen([samtools, "view", "-F", "2308", str(bam)], stdout=subprocess.PIPE, text=True)
    assert p.stdout is not None
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        chrom = f[2]
        targets = marker_positions.get(chrom)
        if not targets:
            continue
        try:
            mapq = int(f[4])
        except ValueError:
            continue
        if mapq < min_mapq:
            continue
        cigar, seq, qual = f[5], f[9].upper(), f[10]
        if cigar == "*" or seq == "*" or qual == "*":
            continue
        for pos, base, bq in base_at_reference_positions(int(f[3]), cigar, seq, qual, targets, min_bq):
            key = (chrom, pos)
            m = lookup.get(key)
            if m is None or base not in {m["REF"], m["ALT"]}:
                continue
            obs[key][base] += 1
            bq_sum[key] += bq
            mapq_sum[key] += mapq
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"samtools view failed while extracting population evidence: {rc}")
    if not obs:
        return _empty_obs(pops)

    rows = []
    for (chrom, pos), counts in sorted(obs.items()):
        m = lookup[(chrom, pos)]
        ref, alt = m["REF"], m["ALT"]
        ref_n, alt_n = int(counts.get(ref, 0)), int(counts.get(alt, 0))
        usable = ref_n + alt_n
        if usable == 0:
            continue
        if ref_n > 0 and alt_n == 0:
            state, allele = "REF", ref
        elif alt_n > 0 and ref_n == 0:
            state, allele = "ALT", alt
        else:
            state, allele = "MIXED", "MIXED"
        row = {
            "Chromosome": chrom, "Position": pos, "REF": ref, "ALT": alt,
            "Observed_state": state, "Observed_allele": allele,
            "REF_observations": ref_n, "ALT_observations": alt_n,
            "Total_usable_base_observations": usable,
            "Mean_BQ": bq_sum[(chrom, pos)] / usable,
            "Mean_MAPQ": mapq_sum[(chrom, pos)] / usable,
            "Block_1Mb": f"{chrom}:{(pos - 1) // block_size_bp}",
            "Max_pairwise_AF_delta": m["Max_pairwise_AF_delta"],
        }
        for pop in pops:
            row[f"AF_{pop}"] = m[f"AF_{pop}"]
        rows.append(row)
    return pd.DataFrame(rows)


def score_population_observations(observations: pd.DataFrame, *, populations=None,
                                  af_floor=0.005, block_size_bp=1_000_000,
                                  bootstraps=1000, seed=20260902,
                                  informative_min_sites=20,
                                  informative_min_bootstrap=0.90) -> dict:
    pops = populations or DEFAULT_POPS
    total = int(len(observations))
    mixed = int((observations["Observed_state"] == "MIXED").sum()) if not observations.empty else 0
    used = observations[observations["Observed_state"].isin(["REF", "ALT"])].copy()
    af_cols = [f"AF_{p}" for p in pops]
    for col in af_cols:
        used[col] = pd.to_numeric(used[col], errors="coerce")
    used = used.dropna(subset=af_cols)
    if used.empty:
        out = {
            "Sites_total_after_collapse": total, "Sites_used": 0,
            "Mixed_sites_excluded": mixed, "Distinct_1Mb_blocks": 0,
            "Best_population": "", "Second_population": "", "Delta_LL": math.nan,
            "Best_bootstrap_support": 0.0, "Population_affinity_recoverable": False,
        }
        for pop in pops:
            out[f"LL_{pop}"] = math.nan
            out[f"Bootstrap_best_{pop}"] = 0.0
        return out

    positions = pd.to_numeric(used["Position"]).astype(int)
    used["Block_1Mb"] = used["Chromosome"].astype(str) + ":" + (((positions - 1) // block_size_bp).astype(str))
    ll_site = np.zeros((len(used), len(pops)), dtype=float)
    states = used["Observed_state"].to_numpy()
    for j, pop in enumerate(pops):
        af = pd.to_numeric(used[f"AF_{pop}"], errors="coerce").to_numpy(dtype=float)
        af = np.clip(af, af_floor, 1.0 - af_floor)
        ll_site[:, j] = np.log(np.where(states == "ALT", af, 1.0 - af))
    ll_total = ll_site.sum(axis=0)
    order = np.argsort(ll_total)[::-1]
    best_i, second_i = int(order[0]), int(order[1]) if len(order) > 1 else int(order[0])

    blocks = used["Block_1Mb"].drop_duplicates().tolist()
    bidx = {b: i for i, b in enumerate(blocks)}
    block_ll = np.zeros((len(blocks), len(pops)), dtype=float)
    for row_i, block in enumerate(used["Block_1Mb"]):
        block_ll[bidx[block], :] += ll_site[row_i, :]
    support = np.zeros(len(pops), dtype=float)
    if bootstraps > 0 and blocks:
        rng = np.random.default_rng(seed)
        counts = np.zeros(len(pops), dtype=int)
        for _ in range(bootstraps):
            draw = rng.integers(0, len(blocks), size=len(blocks))
            counts[int(np.argmax(block_ll[draw, :].sum(axis=0)))] += 1
        support = counts / float(bootstraps)

    best_support = float(support[best_i])
    out = {
        "Sites_total_after_collapse": total,
        "Sites_used": int(len(used)),
        "Mixed_sites_excluded": mixed,
        "Distinct_1Mb_blocks": int(len(blocks)),
        "Best_population": pops[best_i],
        "Second_population": pops[second_i],
        "Delta_LL": float(ll_total[best_i] - ll_total[second_i]),
        "Best_bootstrap_support": best_support,
        "Population_affinity_recoverable": bool(len(used) >= informative_min_sites and best_support >= informative_min_bootstrap),
    }
    for i, pop in enumerate(pops):
        out[f"LL_{pop}"] = float(ll_total[i])
        out[f"Bootstrap_best_{pop}"] = float(support[i])
    return out


def run_population(*, output: str | Path, sample: str, config_path: str | None = None) -> dict:
    output = Path(output).resolve()
    cfg, resolved_cfg = load_effective_config(config_path)
    s = _settings(cfg)
    if not s["enabled"]:
        return {"enabled": False}
    if s["extraction_mode"] != "original_alignment":
        raise RuntimeError("Canonical population extraction_mode must be 'original_alignment'.")

    marker_value = _cfg_get(cfg, "resources", "population_markers", default=None)
    if not marker_value:
        raise RuntimeError("resources.population_markers is not configured")
    marker_path = Path(marker_value)
    if not marker_path.is_file():
        raise FileNotFoundError(f"Population marker resource missing: {marker_path}")

    strict_bam = output / "evidence" / "strict_human.bam"
    if not strict_bam.is_file():
        raise FileNotFoundError(f"Retained strict-human BAM not found: {strict_bam}")

    runtime_value = _cfg_get(cfg, "runtime", "tempdir", default=None)
    runtime_root = Path(runtime_value) if runtime_value else Path(tempfile.gettempdir()) / "srafeguard"
    cache = _marker_cache(marker_path, runtime_root / "population_cache", s["min_pairwise_af_delta"], s["populations"])
    markers = pd.read_csv(cache, sep="\t", compression="gzip")
    samtools = _tool(cfg, "samtools", "samtools")

    obs = _extract_from_original_bam(
        bam=strict_bam, markers=markers, samtools=samtools,
        min_mapq=s["min_mapq"], min_bq=s["min_bq"],
        block_size_bp=s["block_size_bp"], pops=s["populations"],
    )
    obs.to_csv(output / "population_marker_observations.tsv.gz", sep="\t", index=False, compression="gzip")

    score = score_population_observations(
        obs, populations=s["populations"], af_floor=s["af_floor"],
        block_size_bp=s["block_size_bp"], bootstraps=s["bootstraps"], seed=s["seed"],
        informative_min_sites=s["informative_min_sites"],
        informative_min_bootstrap=s["informative_min_bootstrap"],
    )
    score["Sample"] = sample
    cols = ["Sample"] + [c for c in score if c != "Sample"]
    pd.DataFrame([score])[cols].to_csv(output / "population_affinity.tsv", sep="\t", index=False)
    (output / "population_analysis_metadata.json").write_text(json.dumps({
        "sample": sample, "config": resolved_cfg,
        "extraction_mode": "original_alignment_manual_cigar",
        "strict_human_bam": str(strict_bam),
        "population_marker_resource": str(marker_path),
        **{k: v for k, v in s.items() if k != "enabled"},
    }, indent=2))

    print("\n" + "=" * 80)
    print("SRAFEGUARD — POPULATION AFFINITY")
    print("=" * 80)
    print(f"Marker sites observed:     {score['Sites_total_after_collapse']:,}")
    print(f"Informative sites used:   {score['Sites_used']:,}")
    print(f"Mixed sites excluded:     {score['Mixed_sites_excluded']:,}")
    print(f"Distinct 1-Mb blocks:     {score['Distinct_1Mb_blocks']:,}")
    print(f"Best broad affinity:      {score.get('Best_population') or 'None'}")
    print(f"Second broad affinity:    {score.get('Second_population') or 'None'}")
    print(f"Bootstrap support:        {score.get('Best_bootstrap_support', 0.0):.3f}")
    print("Affinity recoverable:     " + ("YES" if score.get("Population_affinity_recoverable") else "NO"))
    return score
