
from pathlib import Path
import tempfile
import yaml

from srafeguard.clinvar import _settings, _resolve_resource_path


def test_hotspot_setting_prefers_canonical_name():
    cfg = {
        "clinvar": {
            "apply_hotspot_exclusion_mask": False,
            "apply_hotspot_mask": True,
        }
    }
    assert _settings(cfg)["apply_hotspot_mask"] is False


def test_hotspot_setting_legacy_alias_still_works():
    cfg = {"clinvar": {"apply_hotspot_mask": False}}
    assert _settings(cfg)["apply_hotspot_mask"] is False


def test_relative_resource_resolution():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "mask.bed"
        p.write_text("chr1\t0\t10\n")
        assert _resolve_resource_path(str(p)) == p.resolve()
