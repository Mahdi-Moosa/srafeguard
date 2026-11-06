from __future__ import annotations
import datetime as dt, functools, hashlib, json, os, shlex, shutil, subprocess, sys, uuid
from pathlib import Path
import pandas as pd, yaml

def _now():
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")

def _repo_root():
    return Path(__file__).resolve().parents[2]

def _plain(v):
    if isinstance(v, Path): return str(v)
    if isinstance(v, dict): return {str(k): _plain(x) for k,x in v.items()}
    if isinstance(v, (list,tuple)): return [_plain(x) for x in v]
    return v

def _atomic_json(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(_plain(obj),indent=2,sort_keys=True,default=str))
    tmp.replace(path)

def _sha256(path,chunk_size=8*1024*1024):
    h=hashlib.sha256()
    with Path(path).open("rb") as fh:
        while True:
            b=fh.read(chunk_size)
            if not b: break
            h.update(b)
    return h.hexdigest()

def _cache_path():
    return Path(os.environ.get("XDG_CACHE_HOME",str(Path.home()/".cache")))/"srafeguard"/"resource_hash_cache.json"

def _load_cache():
    p=_cache_path()
    if not p.is_file(): return {}
    try: return json.loads(p.read_text())
    except Exception: return {}

def _save_cache(x): _atomic_json(_cache_path(),x)

def _resolve_resource(value,config_path=None):
    if not isinstance(value,(str,Path)): return None
    p=Path(value).expanduser()
    if p.is_absolute(): return p
    candidates=[]
    if config_path: candidates.append(Path(config_path).resolve().parent/p)
    candidates += [Path.cwd()/p,_repo_root()/p]
    for c in candidates:
        if c.exists(): return c.resolve()
    return candidates[-1].resolve()

def resource_hashes(cfg,config_path=None):
    resources=cfg.get("resources",{}) or {}
    cache=_load_cache(); changed=False; out={}
    include=("human_reference","population","marker","clinvar_vcf","hotspot","mask")
    exclude=("mmi","index","kraken","database","_db")
    for key,value in resources.items():
        lk=str(key).lower()
        if not any(t in lk for t in include) or any(t in lk for t in exclude):
            continue
        p=_resolve_resource(value,config_path)
        if p is None or not p.is_file():
            out[str(key)]={"path":str(p) if p else str(value),"exists":False,"sha256":None}
            continue
        st=p.stat(); fp=f"{p}|{st.st_size}|{st.st_mtime_ns}"
        digest=cache.get(fp)
        if not digest:
            print(f"+ provenance SHA256: {p}")
            digest=_sha256(p); cache[fp]=digest; changed=True
        out[str(key)]={"path":str(p),"exists":True,"size_bytes":st.st_size,"sha256":digest}
    if changed: _save_cache(cache)
    return out

@functools.lru_cache(maxsize=64)
def _tool_version(executable):
    if not executable: return None
    exe=str(executable)
    if "/" not in exe:
        found=shutil.which(exe)
        if not found: return {"path":None,"version":None}
        exe=found
    elif not Path(exe).exists():
        return {"path":exe,"version":None}
    for suffix in (("--version",),("-version",),("version",)):
        try:
            p=subprocess.run([exe,*suffix],capture_output=True,text=True,timeout=8)
            txt=(p.stdout or p.stderr or "").strip()
            if p.returncode==0 and txt: return {"path":exe,"version":txt.splitlines()[0][:500]}
        except Exception: pass
    return {"path":exe,"version":None}

def tool_versions(cfg):
    return {str(k):_tool_version(v) for k,v in (cfg.get("tools",{}) or {}).items() if isinstance(v,str)}

def _abs(v):
    return str(Path(v).expanduser().resolve()) if v else None

def start_run_record(*,output,sample,cfg,config_path,inputs,mode="audit",parent_run_id=None,context=None):
    output=Path(output).resolve(); output.mkdir(parents=True,exist_ok=True)
    run_id=dt.datetime.now().strftime("%Y%m%d_%H%M%S")+"_"+uuid.uuid4().hex[:8]
    plain_cfg=_plain(cfg)
    effective=output/"effective_config.yaml"
    effective.write_text(yaml.safe_dump(plain_cfg,sort_keys=False,default_flow_style=False,width=120))
    ctx=dict(context or {})
    command=ctx.get("command_invocation") or " ".join(shlex.quote(x) for x in sys.argv)
    rec={
        "run_id":run_id,"parent_run_id":parent_run_id,"mode":mode,"status":"running",
        "sample":sample,"output":str(output),"started_at":_now(),"completed_at":None,
        "command_invocation":command,"cli_overrides":ctx.get("cli_overrides",{}),
        "batch":ctx.get("batch"),"inputs":{k:(_abs(v) if k in {"sra","fastq1","fastq2","reject_source_sra"} else v) for k,v in (inputs or {}).items()},
        "resolved_config_source":str(Path(config_path).resolve()) if config_path else None,
        "effective_config_file":str(effective),"effective_parameters":plain_cfg,
        "tool_versions":tool_versions(cfg),"resource_hashes":resource_hashes(cfg,config_path),"error":None,
    }
    path=output/"run_metadata.json"; _atomic_json(path,rec); return path

def finish_run_record(path,extra=None):
    if not path or not Path(path).is_file(): return
    path=Path(path); rec=json.loads(path.read_text()); rec["status"]="completed"; rec["completed_at"]=_now()
    if extra: rec.update(_plain(extra))
    _atomic_json(path,rec)

def fail_run_record(path,exc,status="failed"):
    if not path or not Path(path).is_file(): return
    path=Path(path)
    try: rec=json.loads(path.read_text())
    except Exception: rec={}
    rec["status"]=status; rec["completed_at"]=_now(); rec["error"]={"type":type(exc).__name__,"message":str(exc)}
    _atomic_json(path,rec)

def annotate_privacy_summary(output,cfg):
    p=Path(output)/"privacy_summary.tsv"
    if not p.is_file(): return
    x=pd.read_csv(p,sep="\t")
    if x.empty: return
    pop=cfg.get("population",{}) or {}; cv=cfg.get("clinvar",{}) or {}
    x["Population_min_MAPQ"]=int(pop.get("min_mapq",30))
    x["Population_min_BQ"]=int(pop.get("min_bq",30))
    x["ClinVar_min_MAPQ"]=int(cv.get("min_mapq",30))
    x["ClinVar_min_BQ"]=int(cv.get("min_bq",30))
    x["ClinVar_chromosome_scope"]=cv.get("chromosome_scope","all")
    x["ClinVar_apply_hotspot_exclusion_mask"]=bool(cv.get("apply_hotspot_exclusion_mask",cv.get("apply_hotspot_mask",True)))
    x.to_csv(p,sep="\t",index=False)
