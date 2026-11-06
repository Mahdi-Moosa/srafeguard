
import pandas as pd

from srafeguard.clinvar import (
    _attach_site_validation,
    _classify_pair_support_rows,
)


def test_pair_validation_rescues_short_alt_mate():
    support = pd.DataFrame([
        {
            "Synthetic_read_id": "x1",
            "QNAME": "q1",
            "Mate": 1,
            "Chromosome": "chr7",
            "Position": 1019,
            "REF": "G",
            "ALT": "A",
        }
    ])

    seq = "C" * 19 + "A" + "C" * 34
    qual = "I" * 54

    pair_maps = {
        "q1": [
            {
                "Mate": 1,
                "Proper_pair": True,
                "Chromosome": "chr7",
                "Start": 1000,
                "End": 1053,
                "MAPQ": 44,
                "CIGAR": "54M",
                "Sequence": seq,
                "Quality": qual,
                "AS": 98,
            },
            {
                "Mate": 2,
                "Proper_pair": True,
                "Chromosome": "chr7",
                "Start": 700,
                "End": 753,
                "MAPQ": 48,
                "CIGAR": "54M",
                "Sequence": "C" * 54,
                "Quality": qual,
                "AS": 108,
            },
        ]
    }

    out = _classify_pair_support_rows(
        support,
        pair_maps,
        min_remap_mapq=30,
        distinct_locus_window_bp=500,
    )

    assert bool(out.iloc[0]["Pair_target_alignment_found"])
    assert bool(out.iloc[0]["Pair_target_ALT_preserved"])
    assert bool(out.iloc[0]["Pair_target_proper_pair"])
    assert bool(out.iloc[0]["Pair_locus_supported"])


def test_site_validation_uses_pair_support_not_independent_only():
    observed = pd.DataFrame([
        {
            "Chromosome": "chr7",
            "Position": 1019,
            "REF": "G",
            "ALT": "A",
            "P_LP_annotation": True,
            "ALT_observations": 1,
        }
    ])

    audit = pd.DataFrame([
        {
            "Chromosome": "chr7",
            "Position": 1019,
            "REF": "G",
            "ALT": "A",
            "QNAME": "q1",
            "Mate": 1,
            "Pair_locus_supported": True,
            "Independent_mate_locus_supported": False,
        }
    ])

    out = _attach_site_validation(
        observed,
        audit,
        min_passing_mates=1,
    )

    assert int(out.iloc[0]["P_LP_pair_locus_supported_mates"]) == 1
    assert int(out.iloc[0]["P_LP_independent_mate_supported_mates"]) == 0
    assert bool(out.iloc[0]["Validated_P_LP_ALT"])
