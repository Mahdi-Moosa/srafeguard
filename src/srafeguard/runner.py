from __future__ import annotations
from pathlib import Path
from .pipeline import run_audit as _pipeline_run_audit
from .provenance import annotate_privacy_summary, fail_run_record, finish_run_record, start_run_record
from .runtime_config import load_effective_config, resolve_config_path

def _infer_sample(*,sample=None,sra=None,fastq1=None):
    if sample: return str(sample)
    if sra:
        n=Path(sra).name
        return n[:-4] if n.lower().endswith(".sra") else n
    if fastq1:
        n=Path(fastq1).name
        for suf in (".fastq.gz",".fq.gz",".fastq",".fq"):
            if n.lower().endswith(suf):
                n=n[:-len(suf)]; break
        for suf in ("_R1","_1",".R1",".1","-R1","-1"):
            if n.endswith(suf):
                n=n[:-len(suf)]; break
        return n or "sample"
    return "sample"

def run_audit(*,sra=None,fastq1=None,fastq2=None,sample=None,reject_source_sra=None,
              layout=None,platform=None,output=None,config_path=None,keep_evidence=None,
              provenance_context=None,**kwargs):
    import os
    from .release_policy import validate_sample,validate_architecture,validate_fastqs,trusted_resources,annotate,clean,diagnostic,SENSITIVE_FILES
    os.umask(0o077)
    resolved_config=resolve_config_path(config_path)
    cfg,_=load_effective_config(resolved_config)
    resolved_sample=_infer_sample(sample=sample,sra=sra,fastq1=fastq1)
    validate_sample(resolved_sample)
    validate_architecture(platform)
    validate_fastqs(fastq1,fastq2)
    trusted_resources(cfg)
    keep_evidence=bool((cfg.get("outputs",{}) or {}).get("keep_evidence",False)) if keep_evidence is None else bool(keep_evidence)
    cfg.setdefault("outputs",{})["keep_evidence"]=keep_evidence
    if output is None:
        output=Path.cwd()/"srafeguard_results"/resolved_sample
    output_path=Path(output).resolve()
    originals=[p for p in (sra,fastq1,fastq2,reject_source_sra) if p]
    owned={output_path/name for name in SENSITIVE_FILES}
    if any(Path(p).resolve() in owned for p in originals):
        raise ValueError('Original input overlaps an owned output; choose a different output directory')
    if output_path.is_symlink():raise ValueError('Output must not be a symlink')
    ctx=dict(provenance_context or {})
    if "cli_overrides" not in ctx:
        vals={"sra":sra,"fastq1":fastq1,"fastq2":fastq2,"sample":sample,
              "reject_source_sra":reject_source_sra,"layout":layout,"platform":platform,
              "output":str(output_path),
              "config":resolved_config,"keep_evidence":keep_evidence}
        ctx["cli_overrides"]={k:v for k,v in vals.items() if v is not None}
    record=None
    analytical_complete=False
    analytical_error=None
    try:
        record=start_run_record(
            output=output_path,sample=resolved_sample,cfg=cfg,config_path=resolved_config,
            inputs={"sra":sra,"fastq1":fastq1,"fastq2":fastq2,
                    "reject_source_sra":reject_source_sra,
                    "layout":layout,"platform":platform},
            mode=ctx.get("mode","audit"),parent_run_id=ctx.get("parent_run_id"),context=ctx)
        result=_pipeline_run_audit(
            sra=sra,fastq1=fastq1,fastq2=fastq2,sample=sample,
            reject_source_sra=reject_source_sra,layout=layout,platform=platform,
            output=str(output_path),
            config_path=resolved_config,keep_evidence=keep_evidence,**kwargs)
        annotate_privacy_summary(output_path,cfg)
        annotate(output_path,cfg)
        analytical_complete=True
        return result
    except KeyboardInterrupt as exc:
        analytical_error=exc
        fail_run_record(record,exc,status="interrupted"); raise
    except BaseException as exc:
        analytical_error=exc
        fail_run_record(record,exc); raise
    finally:
        if record is not None:
            import json
            from .provenance import _atomic_json
            try:
                cleanup=clean(output_path,keep_evidence,protected_paths=originals)
            except Exception as exc:
                cleanup=dict(cleanup_status='FAILED',cleanup_errors=[dict(error_type=type(exc).__name__)],evidence_retention_status='UNKNOWN')
            if analytical_error is not None or cleanup['cleanup_status']!='COMPLETE':
                try:diagnostic(output_path,'analysis' if analytical_error is not None else 'cleanup',type(analytical_error).__name__ if analytical_error is not None else 'CleanupError')
                except OSError:pass
            details=dict(analytical_status='COMPLETED' if analytical_complete else 'INTERRUPTED' if isinstance(analytical_error,KeyboardInterrupt) else 'FAILED',
                         cleanup_status=cleanup['cleanup_status'],cleanup_errors=cleanup['cleanup_errors'],
                         evidence_retention_status=cleanup['evidence_retention_status'],evidence_retention_requested=keep_evidence,
                         evidence_retained=keep_evidence and cleanup['cleanup_status']=='COMPLETE',
                         privacy_summary=str(output_path/'privacy_summary.tsv'))
            if analytical_complete and cleanup['cleanup_status']=='COMPLETE':
                finish_run_record(record,extra=details)
            else:
                d=json.loads(Path(record).read_text());d.update(details)
                if analytical_error is not None and not keep_evidence:
                    d['error']={'type':type(analytical_error).__name__,'message':'Analysis failed; consult the sanitized diagnostic_summary.json.'}
                if analytical_complete:d['status']='cleanup_failed'
                _atomic_json(record,d)
            if analytical_complete and cleanup['cleanup_status']!='COMPLETE':
                raise RuntimeError('Sensitive-output cleanup failed; inspect cleanup status before sharing results')
