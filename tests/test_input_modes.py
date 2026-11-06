from srafeguard.input_modes import (
    canonical_assay,
    canonical_platform,
    infer_layout,
    resolve_input_mode,
)


def test_aliases():
    assert canonical_platform("MGI") == "bgi"
    assert canonical_platform("ONT") == "nanopore"
    assert canonical_assay("RNA-seq") == "transcriptomic"


def test_layout():
    assert infer_layout(fastq1="R1.fq", fastq2="R2.fq") == "paired"
    assert infer_layout(fastq1="reads.fq") == "single"


def test_resolve_short_read():
    cfg = {
        "input": {"assay": "genomic", "platform": "illumina", "layout": "auto"},
        "platforms": {
            "illumina": {
                "enabled": True,
                "read_type": "short",
                "evidence_profile": "short_read",
                "validation_status": "validated_current",
            }
        },
        "assays": {"genomic": {"enabled": True, "splice_aware": False}},
    }
    mode = resolve_input_mode(cfg, fastq1="R1", fastq2="R2")
    assert mode.platform == "illumina"
    assert mode.layout == "paired"
    assert mode.evidence_profile == "short_read"
