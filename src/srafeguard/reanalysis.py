from __future__ import annotations
import copy, datetime as dt, json
from pathlib import Path
import yaml
from .runner import run_audit

def _set_dotted(cfg,dotted,value):
    keys=dotted.split(".")
    if any(not k for k in keys): raise ValueError(f"Invalid dotted config key: {dotted}")
    cur=cfg
    for key in keys[:-1]:
        node=cur.get(key)
        if not isinstance(node,dict):
            node={}; cur[key]=node
        cur=node
    cur[keys[-1]]=value

def _apply_common_overrides(cfg,*,min_mapq=None,min_bq=None,sets=None):
    out=copy.deepcopy(cfg); overrides={}
    if min_mapq is not None:
        for section in ("population","clinvar"): out.setdefault(section,{})["min_mapq"]=int(min_mapq)
        overrides["min_mapq"]=int(min_mapq)
    if min_bq is not None:
        for section in ("population","clinvar"): out.setdefault(section,{})["min_bq"]=int(min_bq)
        overrides["min_bq"]=int(min_bq)
    for item in sets or []:
        if "=" not in item: raise ValueError("--set requires DOTTED.KEY=VALUE")
        key,raw=item.split("=",1); value=yaml.safe_load(raw); _set_dotted(out,key,value); overrides[key]=value
    return out,overrides

def run_reanalysis(*,results,output=None,min_mapq=None,min_bq=None,sets=None,keep_evidence=None,config_path=None):
    parent=Path(results).resolve(); metadata_path=parent/"run_metadata.json"; inherited_config=parent/"effective_config.yaml"
    if not metadata_path.is_file(): raise FileNotFoundError(f"Missing parent run metadata: {metadata_path}")
    if not inherited_config.is_file(): raise FileNotFoundError(f"Missing parent effective config: {inherited_config}")
    metadata=json.loads(metadata_path.read_text()); base_cfg=yaml.safe_load(inherited_config.read_text()) or {}
    if config_path is not None:
        from .release_config import merge
        current=Path(config_path).expanduser().resolve()
        overlay=yaml.safe_load(current.read_text()) or {}
        if not isinstance(overlay,dict):raise ValueError("Current configuration must be a mapping")
        for key,value in overlay.get("resources",{}).items():
            if isinstance(value,str) and value and value!="@bundled_mask":
                p=Path(value).expanduser();overlay["resources"][key]=str(p if p.is_absolute() else (current.parent/p).resolve())
        runtime=overlay.get("runtime",{})
        for key in ("temp_dir","tempdir"):
            if runtime.get(key):
                p=Path(runtime[key]).expanduser();runtime[key]=str(p if p.is_absolute() else (current.parent/p).resolve())
        if "temp_dir" in runtime and "tempdir" not in runtime:runtime["tempdir"]=runtime["temp_dir"]
        if "tempdir" in runtime and "temp_dir" not in runtime:runtime["temp_dir"]=runtime["tempdir"]
        merge(base_cfg,overlay)
    cfg,overrides=_apply_common_overrides(base_cfg,min_mapq=min_mapq,min_bq=min_bq,sets=sets)
    if keep_evidence is None:
        keep_evidence=bool((cfg.get("outputs",{}) or {}).get("keep_evidence",False))
    else:overrides["keep_evidence"]=bool(keep_evidence)
    cfg.setdefault("outputs",{})["keep_evidence"]=bool(keep_evidence)
    child=Path(output).resolve() if output else parent/("reanalyze_"+dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
    child.mkdir(parents=True,exist_ok=False)
    child_cfg=child/"reanalyze_config.yaml"
    child_cfg.write_text(yaml.safe_dump(cfg,sort_keys=False,default_flow_style=False,width=120))
    inputs=metadata.get("inputs") or {}; sample=metadata.get("sample") or parent.name
    return run_audit(
        sra=inputs.get("sra"),fastq1=inputs.get("fastq1"),fastq2=inputs.get("fastq2"),
        sample=sample,reject_source_sra=inputs.get("reject_source_sra"),
        layout=inputs.get("layout"),platform=inputs.get("platform"),
        output=str(child),config_path=str(child_cfg),keep_evidence=keep_evidence,
        provenance_context={"mode":"reanalyze","parent_run_id":metadata.get("run_id"),"cli_overrides":overrides})
