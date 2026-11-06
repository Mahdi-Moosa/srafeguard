from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

import yaml


def _deep_merge(base: dict, override: Mapping[str, Any]) -> dict:
    for key, value in (override or {}).items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def repository_root() -> Path:
    return Path(__file__).resolve().parent

def default_config_path() -> Path:
    return repository_root()/'data'/'default.yaml'

def resolve_config_path(explicit=None):
    raw=explicit or os.environ.get('SRAFEGUARD_CONFIG')
    if not raw:return None
    p=Path(raw).expanduser().resolve()
    if not p.is_file():raise FileNotFoundError(p)
    return str(p)

def load_effective_config(explicit=None):
    from .release_config import load
    resolved=resolve_config_path(explicit)
    return load(resolved),resolved
