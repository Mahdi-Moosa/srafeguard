from pathlib import Path

from srafeguard.clinvar import (
    _is_conservative_plp,
    _merge_bed,
    _reference_end,
    _subtract_bed,
)


def test_reference_end():
    assert _reference_end(100, "100M") == 199
    assert _reference_end(100, "50M2D50M") == 201
    assert _reference_end(100, "5S95M") == 194


def test_merge_and_subtract_bed():
    merged = _merge_bed(
        [
            ("chr1", 10, 20),
            ("chr1", 18, 30),
            ("chr1", 40, 50),
        ]
    )
    assert merged == [
        ("chr1", 10, 30),
        ("chr1", 40, 50),
    ]

    result = _subtract_bed(
        [("chr1", 10, 30)],
        [("chr1", 15, 20)],
    )
    assert result == [
        ("chr1", 10, 15),
        ("chr1", 20, 30),
    ]


def test_conservative_plp():
    assert _is_conservative_plp("Pathogenic")
    assert _is_conservative_plp("Likely_pathogenic")
    assert _is_conservative_plp("Pathogenic/Likely_pathogenic")
    assert not _is_conservative_plp(
        "Conflicting_interpretations_of_pathogenicity"
    )
    assert not _is_conservative_plp("Uncertain_significance")
    assert not _is_conservative_plp("Benign")
