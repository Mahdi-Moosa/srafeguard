"""Fabricated counters only; no real alleles, records or research paths."""
import json
from pathlib import Path
import pandas as pd
import pytest
from srafeguard.release_policy import annotate, clean, SENSITIVE_FILES, SCHEMA
from srafeguard.runtime_config import load_effective_config

def setup_report(tmp_path,layout='PAIRED',**overrides):
    from srafeguard.pipeline import _write_privacy_summary
    from srafeguard.provenance import annotate_privacy_summary
    cfg,_=load_effective_config()
    pd.DataFrame([dict(layout=layout,platform='ARTIFICIAL',evidence_unit='physical_pair' if layout=='PAIRED' else 'single_read',strict_human_fragments=32,strict_unique_qnames=64,chromosomes_with_strict_evidence=2)]).to_csv(tmp_path/'contamination_summary.tsv',sep='\t',index=False)
    cv=dict(ClinVar_validated_P_LP_ALT_sites=2,ClinVar_validated_P_LP_supporting_mates=7,P_LP_validated_exposure_detected=True,ClinVar_sites_with_HQ_base_observation=9,ClinVar_ALT_sites_observed=5,ClinVar_candidate_P_LP_ALT_sites=3)
    if layout=='SINGLE':cv.update(ClinVar_validated_P_LP_supporting_reads=7,ClinVar_validation_mode='single_read_unambiguous_remap',ClinVar_evaluation_status='COMPLETE_SINGLE_READ')
    _write_privacy_summary(output=tmp_path,sample='ARTIFICIAL',input_label='ARTIFICIAL.fastq',population=dict(Sites_used=32,Distinct_1Mb_blocks=4,Best_population='AMR',Second_population='EUR',Best_bootstrap_support=.98,Population_affinity_recoverable=True),clinvar=cv)
    annotate_privacy_summary(tmp_path,cfg)
    old=pd.read_csv(tmp_path/'privacy_summary.tsv',sep='\t',keep_default_na=False)
    for section,value in overrides.items():cfg[section].update(value)
    return cfg,old

def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False)

@pytest.mark.parametrize('layout',['PAIRED','SINGLE'])
def test_exact_legacy_35_fields_and_nonfabricated_counters(tmp_path,layout):
    cfg,old=setup_report(tmp_path,layout);cfg['reporting']['schema']='legacy-1';annotate(tmp_path,cfg)
    legacy=read(tmp_path/'privacy_summary.tsv');pd.testing.assert_frame_equal(old,legacy)
    assert len(legacy.columns)==35 and legacy.ClinVar_validated_P_LP_supporting_mates.iloc[0]==7
    pub=read(tmp_path/'privacy_summary_public.tsv')
    assert pub.Report_schema.iloc[0]==SCHEMA
    if layout=='PAIRED':
        assert pub.Pair_supported_P_LP_ALT_sites.iloc[0]==2 and pub.Pair_supported_P_LP_supporting_mates.iloc[0]==7
        assert pub.Single_read_placement_supporting_reads.iloc[0]=='N/A'
    else:
        assert pub.Pair_supported_P_LP_ALT_sites.iloc[0]=='N/A' and pub.Pair_supported_P_LP_supporting_mates.iloc[0]=='N/A'
        assert pub.Single_read_placement_supporting_reads.iloc[0]==7
        assert legacy.ClinVar_validated_P_LP_supporting_reads.iloc[0]==7

def test_opt_in_sidecar_and_default_public(tmp_path):
    cfg,old=setup_report(tmp_path);cfg['reporting']['legacy_export']=True;annotate(tmp_path,cfg)
    pd.testing.assert_frame_equal(old,read(tmp_path/'privacy_summary_legacy.tsv'))
    assert read(tmp_path/'privacy_summary.tsv').Report_schema.iloc[0]==SCHEMA

def test_no_implicit_legacy_export(tmp_path):
    cfg,_=setup_report(tmp_path);annotate(tmp_path,cfg)
    assert not (tmp_path/'privacy_summary_legacy.tsv').exists()

def test_disabled_modules_keep_configuration_and_missing_not_negative(tmp_path):
    cfg,_=setup_report(tmp_path,clinvar=dict(enabled=False),population=dict(enabled=False));annotate(tmp_path,cfg)
    x=read(tmp_path/'privacy_summary.tsv').iloc[0]
    assert x.ClinVar_ALT_sites_observed=='N/A' and x.Population_affinity_recoverable=='N/A'
    assert x.Population_marker_eligible=='N/A' and not x.Pair_support_evaluated
    for col in ('Population_min_MAPQ','Population_min_BQ','ClinVar_min_MAPQ','ClinVar_min_BQ'):assert x[col]==30
    assert x.ClinVar_chromosome_scope=='all' and isinstance(x.ClinVar_apply_hotspot_exclusion_mask,(bool,__import__('numpy').bool_))
    parsed=pd.read_csv(tmp_path/'privacy_summary.tsv',sep='\t',keep_default_na=False,na_values=['N/A'])
    assert pd.isna(parsed.Population_marker_eligible.iloc[0])
    assert pd.isna(parsed.ClinVar_ALT_sites_observed.astype('Int64').iloc[0])

def test_validation_disabled_is_not_zero(tmp_path):
    cfg,_=setup_report(tmp_path);cfg['clinvar']['validation']['enabled']=False;annotate(tmp_path,cfg)
    x=read(tmp_path/'privacy_summary.tsv').iloc[0]
    assert x.Pair_supported_P_LP_ALT_sites=='N/A' and not x.Pair_support_evaluated
    assert x.ClinVar_ALT_sites_observed==5

@pytest.mark.parametrize('keep',[False,True])
def test_allowlist_aggregates_unknown_files_and_retention(tmp_path,keep):
    for name in SENSITIVE_FILES:
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('ARTIFICIAL')
    aggregates=['population_affinity.tsv','clinvar_exposure.tsv','single_read_clinvar_validation.json','evidence/chromosome_counts.tsv','logs/user_note.txt','evidence/user_note.txt']
    for name in aggregates:(tmp_path/name).write_text('AGGREGATE_OR_USER_FILE')
    result=clean(tmp_path,keep)
    assert result['cleanup_status']=='COMPLETE'
    assert all((tmp_path/name).exists()==keep for name in SENSITIVE_FILES)
    assert all((tmp_path/name).exists() for name in aggregates)

def test_cleanup_errors_explicit_and_original_input_protected(tmp_path,monkeypatch):
    p=tmp_path/'evidence/strict_human.bam';p.parent.mkdir();p.write_text('ORIGINAL')
    result=clean(tmp_path,False,protected_paths=[p]);assert result['cleanup_status']=='FAILED' and p.exists()
    unlink=Path.unlink
    def fail(self,*args,**kwargs):
        if self==p:raise PermissionError('MUST_NOT_EXPORT_RAW_TEXT')
        return unlink(self,*args,**kwargs)
    monkeypatch.setattr(Path,'unlink',fail)
    result=clean(tmp_path,False);assert result['cleanup_errors'][0]['error_type']=='PermissionError'
    assert 'MUST_NOT_EXPORT' not in (tmp_path/'retention_policy.json').read_text()

def test_symlink_directory_not_traversed(tmp_path):
    external=tmp_path/'external';external.mkdir();f=external/'strict_human.bam';f.write_text('ORIGINAL')
    out=tmp_path/'out';out.mkdir();(out/'evidence').symlink_to(external,target_is_directory=True)
    assert clean(out,False)['cleanup_status']=='FAILED' and f.exists()

def fake_run(monkeypatch,tmp_path,fail=False,cleanup_fail=False):
    from srafeguard import runner,provenance,release_policy
    cfg,_=load_effective_config();cfg['release']['reference_profile']='CUSTOM_SYNTHETIC'
    monkeypatch.setattr(runner,'load_effective_config',lambda _: (cfg,None))
    monkeypatch.setattr(provenance,'tool_versions',lambda _: {})
    monkeypatch.setattr(provenance,'resource_hashes',lambda *args: {})
    def pipeline(**kwargs):
        out=Path(kwargs['output']);(out/'logs').mkdir();(out/'logs/human_minimap2.log').write_text('RAW_PRIVATE_ERROR')
        (out/'logs/user.txt').write_text('USER')
        if fail:raise ValueError('PRIVATE_COMMAND_AND_QNAME')
        setup_report(out)
        return {'ok':True}
    monkeypatch.setattr(runner,'_pipeline_run_audit',pipeline)
    if cleanup_fail:
        monkeypatch.setattr(release_policy,'clean',lambda *args,**kwargs:dict(cleanup_status='FAILED',cleanup_errors=[dict(error_type='PermissionError')],evidence_retention_status='PARTIAL_CLEANUP_FAILURE'))
    return runner

def test_runner_failure_sanitized_diagnostic_and_metadata(tmp_path,monkeypatch):
    runner=fake_run(monkeypatch,tmp_path,fail=True)
    with pytest.raises(ValueError):runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out')
    out=tmp_path/'out';d=json.loads((out/'run_metadata.json').read_text())
    assert d['analytical_status']=='FAILED' and d['cleanup_status']=='COMPLETE'
    assert 'PRIVATE_COMMAND' not in (out/'run_metadata.json').read_text()
    assert not (out/'logs/human_minimap2.log').exists() and (out/'logs/user.txt').exists()
    diag=(out/'diagnostic_summary.json').read_text();assert 'ValueError' in diag and 'PRIVATE_COMMAND' not in diag

def test_runner_cleanup_failure_not_marked_complete(tmp_path,monkeypatch):
    runner=fake_run(monkeypatch,tmp_path,cleanup_fail=True)
    with pytest.raises(RuntimeError,match='cleanup failed'):runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out')
    from srafeguard.batch import _completed
    d=json.loads((tmp_path/'out/run_metadata.json').read_text())
    assert d['analytical_status']=='COMPLETED' and d['status']=='cleanup_failed' and d['cleanup_status']=='FAILED'
    assert not _completed(tmp_path/'out')

def test_success_resume_after_cleanup_and_keep_cannot_restore(tmp_path,monkeypatch):
    runner=fake_run(monkeypatch,tmp_path);runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out')
    from srafeguard.batch import _completed
    assert _completed(tmp_path/'out') and not (tmp_path/'out/logs/human_minimap2.log').exists()
    d=json.loads((tmp_path/'out/run_metadata.json').read_text());assert not d['evidence_retained'] and d['cleanup_status']=='COMPLETE'

def test_reanalysis_uses_original_inputs_and_explicit_keep(tmp_path,monkeypatch):
    from srafeguard import reanalysis
    import yaml
    parent=tmp_path/'old';parent.mkdir();original=tmp_path/'original.fastq';original.write_text('ARTIFICIAL')
    (parent/'run_metadata.json').write_text(json.dumps(dict(sample='ARTIFICIAL',inputs=dict(fastq1=str(original),layout='single',platform='ION_TORRENT'))))
    (parent/'effective_config.yaml').write_text(yaml.safe_dump(dict(outputs=dict(keep_evidence=False))))
    seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    reanalysis.run_reanalysis(results=parent,output=tmp_path/'new',keep_evidence=True)
    assert seen['fastq1']==str(original) and seen['keep_evidence'] and original.exists()

def test_invalid_schema_rejected(tmp_path):
    p=tmp_path/'config.yaml';p.write_text('reporting: {schema: unknown}')
    with pytest.raises(ValueError,match='reporting.schema'):load_effective_config(str(p))


def test_analysis_failure_remains_primary_when_cleanup_also_fails(tmp_path,monkeypatch):
    runner=fake_run(monkeypatch,tmp_path,fail=True,cleanup_fail=True)
    with pytest.raises(ValueError,match='PRIVATE_COMMAND'):runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out')
    d=json.loads((tmp_path/'out/run_metadata.json').read_text())
    assert d['analytical_status']=='FAILED' and d['cleanup_status']=='FAILED'
    assert d['status']=='failed' and 'PRIVATE_COMMAND' not in json.dumps(d)

def test_keep_request_preserves_failure_evidence_explicitly(tmp_path,monkeypatch):
    runner=fake_run(monkeypatch,tmp_path,fail=True)
    with pytest.raises(ValueError):runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out',keep_evidence=True)
    d=json.loads((tmp_path/'out/run_metadata.json').read_text())
    assert d['evidence_retention_requested'] and d['evidence_retention_status']=='RETAINED_BY_REQUEST'
    assert (tmp_path/'out/logs/human_minimap2.log').exists()
