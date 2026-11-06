import pandas as pd
from srafeguard.population import base_at_reference_positions, score_population_observations


def test_manual_cigar_projection():
    seq = "TTACGTGCAAA"
    qual = "I" * len(seq)
    hits = base_at_reference_positions(
        100, "2S4M1I3M2D2M", seq, qual,
        [100, 101, 103, 104, 106, 107, 108, 109, 110], 30
    )
    positions = [x[0] for x in hits]
    assert 100 in positions
    assert 103 in positions
    assert 107 not in positions
    assert 108 not in positions


def test_five_population_scoring():
    rows = []
    for i in range(25):
        rows.append({
            "Chromosome":"chr1", "Position":1+i*1100000,
            "REF":"A", "ALT":"G", "Observed_state":"ALT",
            "AF_AFR":0.05, "AF_EUR":0.95, "AF_SAS":0.10,
            "AF_EAS":0.02, "AF_AMR":0.20,
        })
    r = score_population_observations(pd.DataFrame(rows), bootstraps=200)
    assert r["Best_population"] == "EUR"
    assert r["Bootstrap_best_EUR"] == 1.0
    assert r["Population_affinity_recoverable"] is True


def test_mixed_excluded():
    row = {
        "Chromosome":"chr1", "Position":100, "REF":"A", "ALT":"G",
        "Observed_state":"MIXED", "AF_AFR":0.1, "AF_EUR":0.2,
        "AF_SAS":0.3, "AF_EAS":0.4, "AF_AMR":0.5,
    }
    r = score_population_observations(pd.DataFrame([row]))
    assert r["Sites_total_after_collapse"] == 1
    assert r["Mixed_sites_excluded"] == 1
    assert r["Sites_used"] == 0
