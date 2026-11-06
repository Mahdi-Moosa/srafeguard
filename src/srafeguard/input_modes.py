"""Input-mode routing for SRAfeguard.

This module separates routing metadata from validated scientific callers.
Registering a platform here does not itself claim full scientific validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional


PLATFORM_ALIASES = {
    "illumina": "illumina",
    "ilmn": "illumina",
    "bgi": "bgi",
    "mgi": "bgi",
    "dnbseq": "bgi",
    "iontorrent": "iontorrent",
    "ion_torrent": "iontorrent",
    "ion-torrent": "iontorrent",
    "ont": "nanopore",
    "nanopore": "nanopore",
    "oxford_nanopore": "nanopore",
    "pacbio": "pacbio",
    "hifi": "pacbio",
    "pacbio_hifi": "pacbio",
}

ASSAY_ALIASES = {
    "genomic": "genomic",
    "genome": "genomic",
    "wgs": "genomic",
    "dna": "genomic",
    "transcriptomic": "transcriptomic",
    "transcriptome": "transcriptomic",
    "rna": "transcriptomic",
    "rna-seq": "transcriptomic",
    "rnaseq": "transcriptomic",
}

LAYOUT_ALIASES = {
    "auto": "auto",
    "paired": "paired",
    "paired-end": "paired",
    "pe": "paired",
    "single": "single",
    "single-end": "single",
    "se": "single",
    "long": "long",
    "long-read": "long",
}


@dataclass(frozen=True)
class InputMode:
    assay: str
    platform: str
    layout: str
    read_type: str
    evidence_profile: str
    enabled: bool
    validation_status: str
    splice_aware: bool


def _canon(value: Optional[str], aliases: Mapping[str, str], label: str) -> str:
    if value is None:
        return "auto"
    key = str(value).strip().lower()
    if key not in aliases:
        allowed = ", ".join(sorted(set(aliases.values())))
        raise ValueError(f"Unknown {label} '{value}'. Canonical values: {allowed}")
    return aliases[key]


def canonical_platform(value: Optional[str]) -> str:
    if value is None or str(value).strip().lower() == "auto":
        return "auto"
    return _canon(value, PLATFORM_ALIASES, "platform")


def canonical_assay(value: Optional[str]) -> str:
    if value is None or str(value).strip().lower() == "auto":
        return "auto"
    return _canon(value, ASSAY_ALIASES, "assay")


def canonical_layout(value: Optional[str]) -> str:
    if value is None:
        return "auto"
    return _canon(value, LAYOUT_ALIASES, "layout")


def infer_layout(*, fastq1=None, fastq2=None, explicit="auto") -> str:
    layout = canonical_layout(explicit)
    if layout != "auto":
        return layout
    if fastq1 and fastq2:
        return "paired"
    if fastq1 and not fastq2:
        return "single"
    return "auto"


def resolve_input_mode(
    config: Mapping[str, Any],
    *,
    platform: Optional[str] = None,
    assay: Optional[str] = None,
    layout: Optional[str] = None,
    fastq1=None,
    fastq2=None,
    require_enabled: bool = True,
) -> InputMode:
    inp = config.get("input", {}) or {}

    p = canonical_platform(platform or inp.get("platform", "auto"))
    a = canonical_assay(assay or inp.get("assay", "genomic"))
    l = infer_layout(
        fastq1=fastq1,
        fastq2=fastq2,
        explicit=layout or inp.get("layout", "auto"),
    )

    if p == "auto":
        return InputMode(
            assay=a,
            platform=p,
            layout=l,
            read_type="unknown",
            evidence_profile="unresolved",
            enabled=True,
            validation_status="requires_metadata_resolution",
            splice_aware=(a == "transcriptomic"),
        )

    platforms = config.get("platforms", {}) or {}
    assays = config.get("assays", {}) or {}

    if p not in platforms:
        raise ValueError(f"Platform '{p}' has no configuration block")
    if a not in assays:
        raise ValueError(f"Assay '{a}' has no configuration block")

    pcfg = platforms[p] or {}
    acfg = assays[a] or {}

    enabled = bool(pcfg.get("enabled", False)) and bool(acfg.get("enabled", False))
    if require_enabled and not enabled:
        raise ValueError(
            f"Requested mode is configured but disabled: assay={a}, platform={p}. "
            "Enable only after the corresponding scientific caller is validated."
        )

    return InputMode(
        assay=a,
        platform=p,
        layout=l,
        read_type=str(pcfg.get("read_type", "unknown")),
        evidence_profile=str(pcfg.get("evidence_profile", "unresolved")),
        enabled=enabled,
        validation_status=str(pcfg.get("validation_status", "unspecified")),
        splice_aware=bool(acfg.get("splice_aware", False)),
    )
