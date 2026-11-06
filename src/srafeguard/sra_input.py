from pathlib import Path
import shutil
import subprocess

POLICY = "EXCLUDE_REJECT_KEEP_PASS_CRITERIA"

PLATFORM_ENUM = {
    "0": "UNDEFINED",
    "1": "LS454",
    "2": "ILLUMINA",
    "3": "ABI_SOLID",
    "4": "COMPLETE_GENOMICS",
    "5": "HELICOS",
    "6": "PACBIO_SMRT",
    "7": "ION_TORRENT",
    "8": "CAPILLARY",
    "9": "OXFORD_NANOPORE",
    "10": "BGISEQ",
}


def _tool(cfg, name):
    command = str(cfg.get("tools", {}).get(name, name.replace("_", "-")))
    if "/" in command:
        p = Path(command)
        if not p.exists():
            raise FileNotFoundError(f"Tool not found: {p}")
        return str(p)
    found = shutil.which(command)
    if not found:
        raise RuntimeError(f"Required tool not found: {command}")
    return found


def _capture(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}\n"
            + p.stderr[-4000:]
        )
    return p.stdout


def _normalize_layout(value):
    if value is None:
        return None
    s = str(value).strip().upper().replace("-", "_")
    if s in {"", "AUTO"}:
        return None
    if s in {"PAIRED", "PAIRED_END", "PE"}:
        return "PAIRED"
    if s in {"SINGLE", "SINGLE_END", "SE"}:
        return "SINGLE"
    raise ValueError(
        f"Unsupported layout override '{value}'. Use auto, paired, or single."
    )


def _normalize_platform(value):
    if value is None:
        return None
    s = str(value).strip().upper().replace("-", "_").replace(" ", "_")
    if s in {"", "AUTO"}:
        return None
    aliases = {
        "IONTORRENT": "ION_TORRENT",
        "ION_TORRENT": "ION_TORRENT",
        "ILLUMINA": "ILLUMINA",
        "BGI": "BGISEQ",
        "BGISEQ": "BGISEQ",
        "DNBSEQ": "DNBSEQ",
        "MGISEQ": "MGISEQ",
        "PACBIO": "PACBIO_SMRT",
        "PACBIO_SMRT": "PACBIO_SMRT",
        "ONT": "OXFORD_NANOPORE",
        "OXFORD_NANOPORE": "OXFORD_NANOPORE",
    }
    return aliases.get(s, s)


def _platform_from_raw(raw):
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    if s in PLATFORM_ENUM:
        return PLATFORM_ENUM[s]
    s = s.upper().replace("SRA_PLATFORM_", "").replace("-", "_").replace(" ", "_")
    return _normalize_platform(s)


def detect_sra_metadata(sra, cfg, layout=None, platform=None):
    sra = Path(sra)
    layout_override = _normalize_layout(layout)
    platform_override = _normalize_platform(platform)

    cmd = [
        _tool(cfg, "vdb_dump"),
        str(sra),
        "-C",
        "READ_TYPE,PLATFORM",
        "-R",
        "1-50",
        "-f",
        "tab",
    ]

    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode == 0:
        rows = [x for x in p.stdout.splitlines() if x.strip()]
    else:
        txt = _capture(
            [
                _tool(cfg, "vdb_dump"),
                str(sra),
                "-C",
                "READ_TYPE",
                "-R",
                "1-50",
                "-f",
                "tab",
            ]
        )
        rows = [x for x in txt.splitlines() if x.strip()]

    biological_counts = []
    biological_indices = []
    platform_raw = None

    for line in rows:
        f = line.rstrip("\n").split("\t")
        read_type = f[0] if f else ""
        read_types = [x.strip() for x in read_type.split(",")]
        indices = [
            i + 1 for i, typ in enumerate(read_types)
            if "BIOLOGICAL" in typ.upper()
        ]
        if indices:
            biological_counts.append(len(indices))
            if not biological_indices:
                biological_indices = indices
        if len(f) > 1 and not platform_raw:
            candidate = f[1].strip()
            if candidate:
                platform_raw = candidate

    detected_layout = None
    if biological_counts:
        n_bio = max(biological_counts)
        if n_bio >= 2:
            detected_layout = "PAIRED"
        elif n_bio == 1:
            detected_layout = "SINGLE"

    if layout_override:
        resolved_layout = layout_override
        layout_source = "override"
    else:
        resolved_layout = detected_layout
        layout_source = "sra_metadata"

    if not resolved_layout:
        raise RuntimeError(
            "Could not determine SRA library layout from READ_TYPE metadata. "
            "Use --layout paired or --layout single."
        )

    if platform_override:
        resolved_platform = platform_override
        platform_source = "override"
    else:
        resolved_platform = _platform_from_raw(platform_raw) or "UNKNOWN"
        platform_source = "sra_metadata" if platform_raw else "unavailable"

    return {
        "layout": resolved_layout,
        "layout_source": layout_source,
        "platform": resolved_platform,
        "platform_source": platform_source,
        "platform_raw": platform_raw or "",
        "biological_read_indices": biological_indices,
    }


def seqkit_num_seqs(path, cfg):
    txt = _capture([_tool(cfg, "seqkit"), "stats", "-T", str(path)])
    lines = [x for x in txt.splitlines() if x.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"Unexpected seqkit stats output for {path}")
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    return int(float(str(row["num_seqs"]).replace(",", "")))


def sra_reject_spots(sra, reject_file, cfg):
    """
    Frozen WGS_CONTA_2026 policy:
      biological PASS      -> retain
      biological CRITERIA  -> retain
      biological REJECT    -> remove complete spot/evidence unit
    Technical/index REJECT status alone does not remove a biological unit.
    """
    sra = Path(sra)
    reject_file = Path(reject_file)
    reject_file.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        _tool(cfg, "vdb_dump"),
        str(sra),
        "-C",
        "SPOT_ID,READ_TYPE,READ_FILTER",
        "-f",
        "tab",
    ]

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    rejects = []
    rows = 0

    for line in p.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue

        rows += 1
        read_types = [x.strip() for x in f[1].split(",")]
        read_filters = [x.strip() for x in f[2].split(",")]

        biological_reject = False
        for i, flt in enumerate(read_filters):
            typ = read_types[i] if i < len(read_types) else ""
            if "BIOLOGICAL" in typ and "SRA_READ_FILTER_REJECT" in flt:
                biological_reject = True
                break

        if biological_reject:
            rejects.append(int(f[0].strip()))

    stderr = p.stderr.read()
    rc = p.wait()

    if rc != 0:
        raise RuntimeError(
            f"vdb-dump READ_FILTER failed for {sra}: " + stderr[-3000:]
        )
    if rows == 0:
        raise RuntimeError(f"vdb-dump returned no READ_FILTER rows for {sra}")

    with reject_file.open("w") as out:
        for spot in rejects:
            out.write(f"{spot}\n")

    return rejects, rows


def _write_patterns(accession, rejects, path):
    path = Path(path)
    with path.open("w") as out:
        for spot in rejects:
            out.write(f"{accession}.{spot}\n")


def filter_reject_pair(
    accession, r1, r2, reject_spots, work, cfg, require_exact_removal=True
):
    r1 = Path(r1)
    r2 = Path(r2)
    work = Path(work)

    raw1 = seqkit_num_seqs(r1, cfg)
    raw2 = seqkit_num_seqs(r2, cfg)

    if raw1 != raw2:
        raise RuntimeError(
            f"Unequal R1/R2 counts before normalization: {raw1}, {raw2}"
        )

    if not reject_spots:
        return {
            "r1": r1,
            "r2": r2,
            "raw_pairs": raw1,
            "normalized_pairs": raw1,
            "reject_spots": 0,
            "rejects_removed_here": 0,
            "policy": POLICY,
        }

    rf = work / "read_filter"
    rf.mkdir(parents=True, exist_ok=True)
    patterns = rf / "reject_qnames.txt"
    _write_patterns(accession, reject_spots, patterns)

    suffix1 = ".fastq.gz" if str(r1).endswith(".gz") else ".fastq"
    suffix2 = ".fastq.gz" if str(r2).endswith(".gz") else ".fastq"

    out1 = rf / f"{accession}.normalized_1{suffix1}"
    out2 = rf / f"{accession}.normalized_2{suffix2}"

    subprocess.run(
        [_tool(cfg, "seqkit"), "grep", "-v", "-f", str(patterns), str(r1), "-o", str(out1)],
        check=True,
    )
    subprocess.run(
        [_tool(cfg, "seqkit"), "grep", "-v", "-f", str(patterns), str(r2), "-o", str(out2)],
        check=True,
    )

    norm1 = seqkit_num_seqs(out1, cfg)
    norm2 = seqkit_num_seqs(out2, cfg)

    if norm1 != norm2:
        raise RuntimeError(
            f"Unequal R1/R2 counts after normalization: {norm1}, {norm2}"
        )

    removed = raw1 - norm1

    if require_exact_removal:
        if removed != len(reject_spots):
            raise RuntimeError(
                f"SRA normalization mismatch for {accession}: "
                f"expected {len(reject_spots)} rejected pairs, removed {removed}"
            )
    elif removed > len(reject_spots):
        raise RuntimeError(
            f"Removed more reads than documented SRA rejects: "
            f"{removed}>{len(reject_spots)}"
        )

    return {
        "r1": out1,
        "r2": out2,
        "raw_pairs": raw1,
        "normalized_pairs": norm1,
        "reject_spots": len(reject_spots),
        "rejects_removed_here": removed,
        "policy": POLICY,
    }


def filter_reject_single(
    accession, reads, reject_spots, work, cfg, require_exact_removal=True
):
    reads = Path(reads)
    work = Path(work)
    raw = seqkit_num_seqs(reads, cfg)

    if not reject_spots:
        return {
            "r1": reads,
            "r2": None,
            "raw_reads": raw,
            "normalized_reads": raw,
            "reject_spots": 0,
            "rejects_removed_here": 0,
            "policy": POLICY,
        }

    rf = work / "read_filter"
    rf.mkdir(parents=True, exist_ok=True)
    patterns = rf / "reject_qnames.txt"
    _write_patterns(accession, reject_spots, patterns)

    suffix = ".fastq.gz" if str(reads).endswith(".gz") else ".fastq"
    out = rf / f"{accession}.normalized_single{suffix}"

    subprocess.run(
        [_tool(cfg, "seqkit"), "grep", "-v", "-f", str(patterns), str(reads), "-o", str(out)],
        check=True,
    )

    norm = seqkit_num_seqs(out, cfg)
    removed = raw - norm

    if require_exact_removal:
        if removed != len(reject_spots):
            raise RuntimeError(
                f"SRA normalization mismatch for {accession}: "
                f"expected {len(reject_spots)} rejected single reads, removed {removed}"
            )
    elif removed > len(reject_spots):
        raise RuntimeError(
            f"Removed more reads than documented SRA rejects: "
            f"{removed}>{len(reject_spots)}"
        )

    return {
        "r1": out,
        "r2": None,
        "raw_reads": raw,
        "normalized_reads": norm,
        "reject_spots": len(reject_spots),
        "rejects_removed_here": removed,
        "policy": POLICY,
    }


def _select_single_fastq(reads_dir, accession, biological_indices):
    reads_dir = Path(reads_dir)

    direct = reads_dir / f"{accession}.fastq"
    if direct.is_file():
        return direct

    for idx in biological_indices or []:
        p = reads_dir / f"{accession}_{idx}.fastq"
        if p.is_file():
            return p

    produced = sorted(reads_dir.glob("*.fastq"))
    if len(produced) == 1:
        return produced[0]

    candidate = reads_dir / f"{accession}_1.fastq"
    if candidate.is_file() and not biological_indices:
        return candidate

    raise RuntimeError(
        f"Could not unambiguously select the single biological FASTQ for "
        f"{accession}; produced={[x.name for x in produced]}, "
        f"biological_read_indices={biological_indices}"
    )


def extract_and_normalize_sra(
    sra, accession, work, cfg, layout=None, platform=None
):
    sra = Path(sra)
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)

    subprocess.run([_tool(cfg, "vdb_validate"), str(sra)], check=True)

    metadata = detect_sra_metadata(
        sra, cfg, layout=layout, platform=platform
    )

    if metadata["layout"]=="PAIRED" and metadata.get("biological_read_indices") != [1,2]:
        raise ValueError("This SRA paired path only validates biological read indices [1,2]; normalize other spot architectures externally and supply FASTQs.")

    rejects, total_spots = sra_reject_spots(
        sra,
        work / "read_filter" / "reject_spot_ids.txt",
        cfg,
    )

    reads = work / "reads"
    tmp = work / "fasterq_tmp"
    reads.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)

    threads = str(
        cfg.get("input", {}).get("sra", {}).get("fasterq_threads", 12)
    )

    subprocess.run(
        [
            _tool(cfg, "fasterq_dump"),
            str(sra),
            "--split-files",
            "--threads",
            threads,
            "--temp",
            str(tmp),
            "--outdir",
            str(reads),
        ],
        check=True,
    )

    if metadata["layout"] == "PAIRED":
        r1 = reads / f"{accession}_1.fastq"
        r2 = reads / f"{accession}_2.fastq"

        if not r1.exists() or not r2.exists():
            produced = sorted(x.name for x in reads.glob("*.fastq"))
            raise RuntimeError(
                f"Expected paired FASTQs for {accession}; produced={produced}"
            )

        result = filter_reject_pair(
            accession, r1, r2, rejects, work, cfg, require_exact_removal=True
        )
    else:
        single = _select_single_fastq(
            reads,
            accession,
            metadata.get("biological_read_indices", []),
        )
        result = filter_reject_single(
            accession, single, rejects, work, cfg, require_exact_removal=True
        )

    result["total_sra_spots"] = total_spots
    result.update(metadata)
    return result
