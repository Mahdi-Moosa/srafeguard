from __future__ import annotations
import datetime as dt, json
from pathlib import Path
import pandas as pd
from .runner import run_audit
from .runtime_config import resolve_config_path

def _clean(v):
    if v is None: return None
    s=str(v).strip(); return s or None

def load_manifest(path):
    df=pd.read_csv(Path(path),sep="\t",dtype=str).fillna("")
    if "sample" not in df.columns: raise ValueError("Batch manifest requires a 'sample' column.")
    if not ("sra" in df.columns or "fastq1" in df.columns):
        raise ValueError("Batch manifest requires an 'sra' or 'fastq1' column.")
    if df["sample"].duplicated().any():
        dup=df.loc[df["sample"].duplicated(),"sample"].tolist()
        raise ValueError(f"Duplicate sample identifiers: {dup[:10]}")
    from .release_policy import validate_sample
    for sample in df["sample"]:validate_sample(sample)
    return df

def _completed(output_dir):
    o=Path(output_dir); meta=o/"run_metadata.json"; summary=o/"privacy_summary.tsv"
    if not (meta.is_file() and summary.is_file()): return False
    try:
        record=json.loads(meta.read_text())
        return record.get("status")=="completed" and record.get("analytical_status")=="COMPLETED" and record.get("cleanup_status")=="COMPLETE"
    except Exception: return False

def run_batch(*,manifest,output_root,config_path=None,keep_evidence=None,resume=True,fail_fast=False):
    manifest=Path(manifest).resolve(); output_root=Path(output_root).resolve(); output_root.mkdir(parents=True,exist_ok=True)
    config_path=resolve_config_path(config_path); df=load_manifest(manifest); n=len(df)
    batch_id=dt.datetime.now().strftime("%Y%m%d_%H%M%S")+"_"+manifest.stem
    rows=[]
    for i,row in df.iterrows():
        sample=_clean(row.get("sample")); out=output_root/sample
        if resume and _completed(out):
            print(f"[{i+1}/{n}] {sample} — SKIP completed")
            rows.append({"sample":sample,"status":"skipped_completed","output":str(out),"error":""})
            continue
        sra=_clean(row.get("sra")); fastq1=_clean(row.get("fastq1")); fastq2=_clean(row.get("fastq2"))
        reject=_clean(row.get("reject_source_sra"))
        layout=_clean(row.get("layout")) or "auto"
        platform=_clean(row.get("platform")) or "auto"
        if bool(sra)==bool(fastq1):
            status="failed"; err="Each manifest row must define exactly one of 'sra' or 'fastq1'."
        elif fastq1 and layout.lower() == "paired" and not fastq2:
            status="failed"; err="'fastq2' is required for paired FASTQ input."
        elif fastq1 and layout.lower() == "single" and fastq2:
            status="failed"; err="'fastq2' cannot be supplied for single FASTQ input."
        else:
            print(); print("="*80); print(f"[{i+1}/{n}] {sample}"); print("="*80)
            try:
                run_audit(
                    sra=sra,fastq1=fastq1,fastq2=fastq2,sample=sample,reject_source_sra=reject,
                    layout=layout,platform=platform,
                    output=str(out),config_path=config_path,keep_evidence=keep_evidence,
                    provenance_context={
                        "mode":"batch_audit",
                        "cli_overrides":{"keep_evidence":keep_evidence},
                        "batch":{"batch_id":batch_id,"manifest":str(manifest),
                                 "row_index_1based":i+1,"total_rows":n}})
                status="completed"; err=""
            except KeyboardInterrupt:
                raise
            except BaseException as exc:
                status="failed"; err=str(exc)
        print(f"[{i+1}/{n}] {sample} — {status.upper()}")
        if err: print(err)
        rows.append({"sample":sample,"status":status,"output":str(out),"error":err})
        pd.DataFrame(rows).to_csv(output_root/"batch_status.tsv",sep="\t",index=False)
        if status=="failed" and fail_fast: break
    result=pd.DataFrame(rows)
    result.to_csv(output_root/"batch_status.tsv",sep="\t",index=False)
    print(); print("Batch status:")
    if not result.empty: print(result["status"].value_counts().to_string())
    return result
