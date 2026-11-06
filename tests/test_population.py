import pandas as pd

from srafeguard.population import score_population_observations


def test_five_population_scoring_recoverable_eur():
    pops = ["AFR", "EUR", "SAS", "EAS", "AMR"]
    rows = []

    for i in range(25):
        row = {
            "Chromosome": "chr1",
            "Position": 1 + i * 1_100_000,
            "REF": "A",
            "ALT": "G",
            "Observed_state": "ALT",
            "AF_AFR": 0.05,
            "AF_EUR": 0.95,
            "AF_SAS": 0.10,
            "AF_EAS": 0.02,
            "AF_AMR": 0.20,
        }
        rows.append(row)

    result = score_population_observations(
        pd.DataFrame(rows),
        populations=pops,
        bootstraps=200,
        seed=20260902,
        informative_min_sites=20,
        informative_min_bootstrap=0.90,
    )

    assert result["Sites_used"] == 25
    assert result["Best_population"] == "EUR"
    assert result["Bootstrap_best_EUR"] == 1.0
    assert result["Population_affinity_recoverable"] is True


def test_population_scoring_empty_is_not_recoverable():
    cols = [
        "Chromosome",
        "Position",
        "REF",
        "ALT",
        "Observed_state",
        "AF_AFR",
        "AF_EUR",
        "AF_SAS",
        "AF_EAS",
        "AF_AMR",
    ]
    result = score_population_observations(
        pd.DataFrame(columns=cols)
    )

    assert result["Sites_used"] == 0
    assert result["Best_population"] == ""
    assert result["Population_affinity_recoverable"] is False
