
import pandas as pd

from srafeguard.clinvar import (
    _classify_plp_support_rows,
    _filter_territory_scope,
)


def test_autosome_scope():
    x = [
        ("chr1", 0, 10),
        ("chrX", 0, 10),
        ("chrY", 0, 10),
        ("chrM", 0, 10),
    ]
    assert _filter_territory_scope(x, "autosomes") == [
        ("chr1", 0, 10)
    ]
    assert _filter_territory_scope(x, "all") == x


def test_generic_plp_validator_accepts_unique_target_alignment():
    support = pd.DataFrame(
        [
            {
                "Synthetic_read_id": "r1",
                "QNAME": "q1",
                "Mate": 1,
                "Chromosome": "chr7",
                "Position": 1005,
                "REF": "G",
                "ALT": "A",
                "ClinVar_ID": "x",
            }
        ]
    )

    maps = {
        "r1": [
            {
                "Chromosome": "chr7",
                "Start": 1000,
                "End": 1099,
                "MAPQ": 60,
                "CIGAR": "100M",
                "Sequence": "AAAAAA",
                "Quality": "IIIIII",
                "AS": 100,
            },
            {
                "Chromosome": "chr12",
                "Start": 5000,
                "End": 5099,
                "MAPQ": 0,
                "CIGAR": "100M",
                "Sequence": "*",
                "Quality": "*",
                "AS": 70,
            },
        ]
    }

    # Replace the sequence with 100 bp and ALT A at target offset 5.
    maps["r1"][0]["Sequence"] = "A" * 100
    maps["r1"][0]["Quality"] = "I" * 100

    out = _classify_plp_support_rows(
        support,
        maps,
        min_remap_mapq=30,
        distinct_locus_window_bp=500,
    )

    assert bool(out.iloc[0]["Target_alignment_found"])
    assert bool(out.iloc[0]["Target_ALT_preserved_on_remap"])
    assert bool(out.iloc[0]["Unambiguous_locus_support"])


def test_generic_plp_validator_rejects_equal_competitor():
    support = pd.DataFrame(
        [
            {
                "Synthetic_read_id": "r1",
                "QNAME": "q1",
                "Mate": 1,
                "Chromosome": "chr7",
                "Position": 1005,
                "REF": "G",
                "ALT": "A",
                "ClinVar_ID": "x",
            }
        ]
    )

    seq = "A" * 100
    qual = "I" * 100

    maps = {
        "r1": [
            {
                "Chromosome": "chr7",
                "Start": 1000,
                "End": 1099,
                "MAPQ": 30,
                "CIGAR": "100M",
                "Sequence": seq,
                "Quality": qual,
                "AS": 100,
            },
            {
                "Chromosome": "chr12",
                "Start": 5000,
                "End": 5099,
                "MAPQ": 0,
                "CIGAR": "100M",
                "Sequence": "*",
                "Quality": "*",
                "AS": 100,
            },
        ]
    }

    out = _classify_plp_support_rows(
        support,
        maps,
        min_remap_mapq=30,
        distinct_locus_window_bp=500,
    )

    assert not bool(out.iloc[0]["Unambiguous_locus_support"])
