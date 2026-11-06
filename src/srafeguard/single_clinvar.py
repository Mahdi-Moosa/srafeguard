from __future__ import annotations

import ast
import contextlib
import inspect
import io
import json
import textwrap
from pathlib import Path
import pandas as pd

from . import clinvar as cv

VALIDATION_MODE = "single_read_unambiguous_remap"


def _extract_merge_keys(fn):
    try:
        src = textwrap.dedent(inspect.getsource(fn))
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                if any(isinstance(t, ast.Name) and t.id == "keys" for t in node.targets):
                    value = ast.literal_eval(node.value)
                    if isinstance(value, list) and all(isinstance(x, str) for x in value):
                        return value
    except Exception:
        pass
    return ["Chromosome", "Position", "ALT", "QNAME", "Mate"]


def _touch_pair_outputs(*args, **kwargs):
    for key in ("out_r1", "out_r2"):
        value = kwargs.get(key)
        if value is not None:
            p = Path(value)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.touch()
    return None


def _empty_pair_maps(*args, **kwargs):
    return {}


def _empty_pair_audit(*args, **kwargs):
    return pd.DataFrame()


@contextlib.contextmanager
def _single_read_validation_context():
    required = [
        "_extract_supporting_pairs_from_bam",
        "_remap_pairs",
        "_classify_pair_support_rows",
        "_merge_independent_and_pair_validation",
        "_classify_plp_support_rows",
    ]
    missing = [name for name in required if not hasattr(cv, name)]
    if missing:
        raise RuntimeError(
            "Installed ClinVar module is incompatible with the single-read adapter; "
            f"missing functions: {missing}"
        )

    original_extract = cv._extract_supporting_pairs_from_bam
    original_remap_pairs = cv._remap_pairs
    original_classify_pair = cv._classify_pair_support_rows
    original_merge = cv._merge_independent_and_pair_validation
    merge_keys = _extract_merge_keys(original_merge)

    def merge_single(independent, pair_validation):
        if independent is None or independent.empty:
            return original_merge(
                independent if independent is not None else pd.DataFrame(),
                pd.DataFrame(),
            )

        if "Unambiguous_locus_support" not in independent.columns:
            raise RuntimeError(
                "Independent ClinVar remap audit lacks Unambiguous_locus_support."
            )

        keys = [x for x in merge_keys if x in independent.columns]
        if not keys:
            raise RuntimeError(
                f"Could not identify ClinVar merge keys from {merge_keys}."
            )

        synthetic = independent[keys].copy()
        synthetic["Pair_locus_supported"] = (
            independent["Unambiguous_locus_support"]
            .fillna(False)
            .astype(bool)
            .to_numpy()
        )
        synthetic["Pair_validation_mode"] = VALIDATION_MODE
        return original_merge(independent, synthetic)

    cv._extract_supporting_pairs_from_bam = _touch_pair_outputs
    cv._remap_pairs = _empty_pair_maps
    cv._classify_pair_support_rows = _empty_pair_audit
    cv._merge_independent_and_pair_validation = merge_single
    try:
        yield
    finally:
        cv._extract_supporting_pairs_from_bam = original_extract
        cv._remap_pairs = original_remap_pairs
        cv._classify_pair_support_rows = original_classify_pair
        cv._merge_independent_and_pair_validation = original_merge


def run_clinvar_single(*, output, sample, config_path=None):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        with _single_read_validation_context():
            summary = cv.run_clinvar(
                output=output,
                sample=sample,
                config_path=config_path,
            )

    text = buf.getvalue().replace(
        "Computational support mates:",
        "Computational support reads:",
    )
    if text:
        print(text, end="")

    summary = dict(summary or {})
    summary["ClinVar_evaluation_status"] = "COMPLETE_SINGLE_READ"
    summary["ClinVar_validation_mode"] = VALIDATION_MODE
    summary["ClinVar_validated_P_LP_supporting_reads"] = summary.get(
        "ClinVar_validated_P_LP_supporting_mates", 0
    )

    record = {
        "validation_mode": VALIDATION_MODE,
        "evidence_unit": "single_read",
        "paired_path_modified": False,
        "canonical_rule": (
            "ALT-supporting read remaps over the target locus with ALT preserved, "
            "meets the installed remap MAPQ threshold, and has no disqualifying "
            "equal-or-better distinct-locus competitor."
        ),
        "validated_plp_alt_sites": summary.get(
            "ClinVar_validated_P_LP_ALT_sites", 0
        ),
        "validated_supporting_reads": summary.get(
            "ClinVar_validated_P_LP_supporting_reads", 0
        ),
    }
    Path(output, "single_read_clinvar_validation.json").write_text(
        json.dumps(record, indent=2) + "\n"
    )
    return summary
