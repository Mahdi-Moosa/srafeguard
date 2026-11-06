from pathlib import Path
import shutil

from .config import (
    PACKAGE_ROOT,
    get_dotted,
    load_config,
    load_yaml,
)


RESOURCE_MANIFEST = PACKAGE_ROOT / "data" / "resources.yaml"


def human_size(n):
    if n is None:
        return "-"

    x = float(n)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024 or unit == "TB":
            return f"{x:.2f} {unit}"
        x /= 1024


def resolve_resource_path(name, spec, cfg):
    if spec.get("type") == "bundled":
        rel = spec.get("path")
        return PACKAGE_ROOT / rel if rel else None

    key = spec.get("config_key")

    if not key:
        return None

    value = get_dotted(cfg, key)

    if not value:
        return None

    p = Path(value)

    if not p.is_absolute():
        p = PACKAGE_ROOT / p

    return p


def check_resources(config_path=None):
    cfg = load_config(config_path)
    manifest = load_yaml(RESOURCE_MANIFEST)

    rows = []

    for name, spec in manifest.get("resources", {}).items():
        path = resolve_resource_path(name, spec, cfg)

        if path is None:
            status = "UNCONFIGURED"
            size = None
        elif not path.exists():
            status = "MISSING"
            size = None
        elif not path.is_file():
            status = "NOT_FILE"
            size = None
        elif path.stat().st_size == 0:
            status = "EMPTY"
            size = 0
        else:
            status = "OK"
            size = path.stat().st_size

        rows.append(
            {
                "resource": name,
                "status": status,
                "size": size,
                "path": str(path) if path else "",
                "required": bool(spec.get("required", False)),
                "type": spec.get("type", ""),
            }
        )

    return cfg, rows


def check_tools(cfg):
    rows = []

    for name, command in cfg.get("tools", {}).items():
        if not command:
            rows.append((name, "UNCONFIGURED", ""))
            continue

        command = str(command)

        if "/" in command:
            p = Path(command)
            found = str(p) if p.exists() else None
        else:
            found = shutil.which(command)

        rows.append(
            (
                name,
                "OK" if found else "MISSING",
                found or command,
            )
        )

    return rows


def print_resource_check(config_path=None):
    cfg, rows = check_resources(config_path)

    print("=" * 118)
    print("SRAFEGUARD — RESOURCE CHECK")
    print("=" * 118)

    for r in rows:
        print(
            f"{r['resource']:28s} "
            f"{r['status']:13s} "
            f"{human_size(r['size']):>11s}  "
            f"{r['path']}"
        )

    print()
    print("TOOLS")
    print("-" * 118)

    for name, status, location in check_tools(cfg):
        print(
            f"{name:28s} "
            f"{status:13s} "
            f"{location}"
        )

    required_bad = [
        r for r in rows
        if r["required"] and r["status"] != "OK"
    ]

    print()
    print("-" * 118)

    if required_bad:
        print(
            f"RESOURCE CHECK FAILED: "
            f"{len(required_bad)} required resource(s) unavailable."
        )
        return 1

    print("RESOURCE CHECK PASSED")
    return 0
