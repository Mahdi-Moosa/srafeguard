from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


def _find_existing_bam(obj: Any) -> Path | None:
    if isinstance(obj, (str, os.PathLike)):
        p = Path(obj)
        if p.suffix.lower() == ".bam" and p.is_file():
            return p.resolve()
    if isinstance(obj, dict):
        for value in obj.values():
            p = _find_existing_bam(value)
            if p is not None:
                return p
    if isinstance(obj, (list, tuple, set)):
        for value in obj:
            p = _find_existing_bam(value)
            if p is not None:
                return p
    return None


def _preserve_path(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.unlink(missing_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


@contextlib.contextmanager
def capture_original_human_bam(destination: str | Path):
    """Retain BAM produced by either validated human mapper."""
    from . import paired_caller_frozen
    from . import single_caller_frozen

    destination = Path(destination).resolve()
    paired_original = paired_caller_frozen.human_map
    single_original = single_caller_frozen.human_map_single

    def _wrap(original):
        def wrapped(*args, **kwargs):
            result = original(*args, **kwargs)
            bam = (
                _find_existing_bam(result)
                or _find_existing_bam(args)
                or _find_existing_bam(kwargs)
            )
            if bam is None:
                raise RuntimeError(
                    "Could not identify the BAM produced by the validated human mapper."
                )
            _preserve_path(bam, destination)
            return result
        return wrapped

    paired_caller_frozen.human_map = _wrap(paired_original)
    single_caller_frozen.human_map_single = _wrap(single_original)

    try:
        yield destination
    finally:
        paired_caller_frozen.human_map = paired_original
        single_caller_frozen.human_map_single = single_original


def finalize_strict_human_bam(
    *,
    captured_bam: str | Path,
    qnames_file: str | Path,
    output_bam: str | Path,
    samtools: str,
    sort_threads: int = 2,
) -> Path:
    captured_bam = Path(captured_bam)
    qnames_file = Path(qnames_file)
    output_bam = Path(output_bam)

    if not captured_bam.is_file():
        raise FileNotFoundError(f"Captured original human BAM not found: {captured_bam}")
    if not qnames_file.is_file():
        raise FileNotFoundError(f"Strict-QNAME file not found: {qnames_file}")

    output_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_bam.with_suffix(".selected.unsorted.bam")

    from .owned_temporary import register
    import tempfile
    output=output_bam.parent.parent
    register(output,tmp,'strict_selection')
    sort_work=Path(tempfile.mkdtemp(prefix='srafeguard_sort_',dir=output_bam.parent))
    register(output,sort_work,'sort_workspace')
    try:
        subprocess.run(
            [samtools, 'view', '-b', '-F', '2308', '-N', str(qnames_file),
             '-o', str(tmp), str(captured_bam)], check=True)
        subprocess.run(
            [samtools, 'sort', '-@', str(sort_threads), '-T', str(sort_work/'chunk'),
             '-o', str(output_bam), str(tmp)], check=True)
        subprocess.run([samtools, 'index', str(output_bam)], check=True)
        return output_bam
    finally:
        # The outer ownership cleanup retries failures and reports any remaining data.
        try:tmp.unlink(missing_ok=True)
        except OSError:pass
        try:shutil.rmtree(sort_work)
        except OSError:pass
