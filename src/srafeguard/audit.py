from pathlib import Path
import shutil
import tempfile

from .config import load_config
from .contamination import run_contamination
from .sra_input import (
    extract_and_normalize_sra,
    filter_reject_pair,
    filter_reject_single,
    sra_reject_spots,
    _normalize_layout,
    _normalize_platform,
)


def infer_sample_from_sra(sra):
    name = Path(sra).name
    if name.endswith(".sra"):
        name = name[:-4]
    return name


def run_audit(
    output,
    config_path=None,
    sample=None,
    sra=None,
    fastq1=None,
    fastq2=None,
    reject_source_sra=None,
    layout=None,
    platform=None,
):
    cfg = load_config(config_path)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    temp_root = cfg.get("runtime", {}).get("temp_dir")
    if temp_root:
        Path(temp_root).mkdir(parents=True, exist_ok=True)

    work = Path(
        tempfile.mkdtemp(prefix="hla_input_", dir=temp_root)
    )

    from .owned_temporary import register
    register(output,work,"input_workspace")

    try:
        if sra:
            if fastq1 or fastq2:
                raise ValueError("--sra cannot be combined with --fastq1/--fastq2")

            sample = sample or infer_sample_from_sra(sra)
            inp = extract_and_normalize_sra(
                sra=sra,
                accession=sample,
                work=work,
                cfg=cfg,
                layout=layout,
                platform=platform,
            )

            r1 = inp["r1"]
            r2 = inp.get("r2")
            resolved_layout = inp["layout"]
            resolved_platform = inp["platform"]

            print(
                f"+ SRA metadata: layout={resolved_layout} ({inp['layout_source']}), "
                f"platform={resolved_platform} ({inp['platform_source']})"
            )

            if resolved_layout == "PAIRED":
                print(
                    f"+ SRA READ_FILTER normalization: "
                    f"{inp['reject_spots']} rejected biological pair(s); "
                    f"{inp['normalized_pairs']:,} pairs retained"
                )
            else:
                print(
                    f"+ SRA READ_FILTER normalization: "
                    f"{inp['reject_spots']} rejected biological read(s); "
                    f"{inp['normalized_reads']:,} reads retained"
                )
        else:
            if not fastq1:
                raise ValueError("Provide either --sra or --fastq1.")

            if not sample:
                sample = Path(fastq1).name.split("_")[0]

            r1 = Path(fastq1)
            r2 = Path(fastq2) if fastq2 else None

            override_layout = _normalize_layout(layout)
            resolved_layout = override_layout or ("PAIRED" if r2 is not None else "SINGLE")

            if resolved_layout == "PAIRED" and r2 is None:
                raise ValueError("PAIRED layout requires --fastq2.")
            if resolved_layout == "SINGLE" and r2 is not None:
                raise ValueError("SINGLE layout cannot be combined with --fastq2.")

            resolved_platform = _normalize_platform(platform) or "UNKNOWN"

            if reject_source_sra:
                rejects, _ = sra_reject_spots(
                    reject_source_sra,
                    work / "read_filter" / "reject_spot_ids.txt",
                    cfg,
                )

                if resolved_layout == "PAIRED":
                    inp = filter_reject_pair(
                        sample, r1, r2, rejects, work, cfg,
                        require_exact_removal=False,
                    )
                    r1 = inp["r1"]
                    r2 = inp["r2"]
                    print(
                        f"+ Post-filter SRA normalization: "
                        f"{inp['rejects_removed_here']} surviving REJECT pair(s) removed; "
                        f"{inp['normalized_pairs']:,} pairs retained"
                    )
                else:
                    inp = filter_reject_single(
                        sample, r1, rejects, work, cfg,
                        require_exact_removal=False,
                    )
                    r1 = inp["r1"]
                    print(
                        f"+ Post-filter SRA normalization: "
                        f"{inp['rejects_removed_here']} surviving REJECT read(s) removed; "
                        f"{inp['normalized_reads']:,} reads retained"
                    )

        return run_contamination(
            fastq1=r1,
            fastq2=r2,
            sample=sample,
            output=output,
            config_path=config_path,
            layout=resolved_layout,
            platform=resolved_platform,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)
