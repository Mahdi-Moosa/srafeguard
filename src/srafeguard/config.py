from copy import deepcopy
from pathlib import Path
import yaml


PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PACKAGE_ROOT / "data" / "default.yaml"


def load_yaml(path):
    path = Path(path)
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def deep_merge(base, override):
    out = deepcopy(base)

    for key, value in override.items():
        if (
            key in out
            and isinstance(out[key], dict)
            and isinstance(value, dict)
        ):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = deepcopy(value)

    return out


def load_config(local_config=None):
    from .runtime_config import load_effective_config
    return load_effective_config(local_config)[0]


def get_dotted(mapping, dotted_key, default=None):
    cur = mapping

    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]

    return cur
