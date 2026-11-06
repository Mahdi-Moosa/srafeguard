"""Fabricated inputs only: exceptional cleanup and retention-precedence gates."""
import json,subprocess,sys,io
from pathlib import Path
from unittest.mock import patch
import pytest,yaml
from srafeguard.alignment_evidence import finalize_strict_human_bam
from srafeguard.release_policy import clean
from srafeguard import reanalysis,cli

def finalize_fixture(p):
    out=p/'out';(out/'evidence').mkdir(parents=True)
    src=p/'ARTIFICIAL.bam';src.write_bytes(b'ARTIFICIAL_INPUT')
    q=p/'ARTIFICIAL.qnames';q.write_text('ARTIFICIAL\n')
    return out,src,q

@pytest.mark.parametrize('failure',['view','sort','exception','interrupt'])
def test_owned_temporary_absent_after_default_failure_cleanup(tmp_path,failure):
    out,src,q=finalize_fixture(tmp_path);dest=out/'evidence/strict_human.bam'
    def execute(args,**kw):
        if args[1]=='view':
            Path(args[args.index('-o')+1]).write_bytes(b'ARTIFICIAL_TEMP')
            if failure=='view':raise subprocess.CalledProcessError(1,args)
        if args[1]=='sort':
            Path(args[args.index('-o')+1]).write_bytes(b'ARTIFICIAL_PARTIAL')
            if failure=='interrupt':raise KeyboardInterrupt()
            if failure=='exception':raise OSError('ARTIFICIAL_SORT_EXCEPTION')
            raise subprocess.CalledProcessError(1,args)
        return subprocess.CompletedProcess(args,0)
    with patch('srafeguard.alignment_evidence.subprocess.run',side_effect=execute):
        with pytest.raises((subprocess.CalledProcessError,OSError,KeyboardInterrupt)):
            finalize_strict_human_bam(captured_bam=src,qnames_file=q,output_bam=dest,samtools='ARTIFICIAL')
    state=clean(out,False,protected_paths=[src,q])
    assert state['cleanup_status']=='COMPLETE'
    assert not dest.exists() and not dest.with_suffix('.selected.unsorted.bam').exists()
    assert src.read_bytes()==b'ARTIFICIAL_INPUT' and q.exists()

def parent(p,keep):
    p.mkdir();(p/'run_metadata.json').write_text(json.dumps(dict(sample='ARTIFICIAL',inputs=dict(fastq1='ARTIFICIAL.fastq',layout='single'))))
    (p/'effective_config.yaml').write_text(yaml.safe_dump(dict(outputs=dict(keep_evidence=keep))))
    return p

@pytest.mark.parametrize('old,new',[(True,False),(False,True)])
def test_reanalysis_set_honored_with_effective_configuration(tmp_path,monkeypatch,old,new):
    p=parent(tmp_path/'parent',old);seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    reanalysis.run_reanalysis(results=p,output=tmp_path/'child',sets=['outputs.keep_evidence='+str(new).lower()])
    child=yaml.safe_load(Path(seen['config_path']).read_text())
    assert seen['keep_evidence']==new and child['outputs']['keep_evidence']==new

@pytest.mark.parametrize('old,new',[(True,False),(False,True)])
def test_reanalysis_explicit_cli_has_highest_precedence(tmp_path,monkeypatch,old,new):
    p=parent(tmp_path/'parent',old);seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    reanalysis.run_reanalysis(results=p,output=tmp_path/'child',sets=['outputs.keep_evidence='+str(old).lower()],keep_evidence=new)
    assert seen['keep_evidence']==new
    assert yaml.safe_load(Path(seen['config_path']).read_text())['outputs']['keep_evidence']==new

@pytest.mark.parametrize('old,new',[(True,False),(False,True)])
def test_current_config_overrides_inherited_retention(tmp_path,monkeypatch,old,new):
    p=parent(tmp_path/'parent',old);seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    config=tmp_path/'current.yaml';config.write_text(yaml.safe_dump(dict(outputs=dict(keep_evidence=new),resources=dict(population_markers='external/ARTIFICIAL.tsv.gz'))))
    reanalysis.run_reanalysis(results=p,output=tmp_path/'child',config_path=config)
    cfg=yaml.safe_load(Path(seen['config_path']).read_text())
    assert seen['keep_evidence']==new and cfg['outputs']['keep_evidence']==new
    assert cfg['resources']['population_markers']==str(tmp_path/'external/ARTIFICIAL.tsv.gz')

@pytest.mark.parametrize('keep',[False,True])
def test_no_override_inherits_without_copying_evidence(tmp_path,monkeypatch,keep):
    p=parent(tmp_path/'parent',keep);(p/'evidence').mkdir();(p/'evidence/private.txt').write_text('ARTIFICIAL_OLD')
    seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    reanalysis.run_reanalysis(results=p,output=tmp_path/'child')
    assert seen['keep_evidence']==keep and not (tmp_path/'child/evidence').exists()
    assert (p/'evidence/private.txt').exists()

def test_all_cli_retention_options_are_tristate():
    parser=cli.build_parser()
    for base in [['audit','--fastq1','ARTIFICIAL.fq'],['batch','--manifest','ARTIFICIAL.tsv','--output-root','ARTIFICIAL'],['reanalyze','--results','ARTIFICIAL']]:
        assert parser.parse_args(base).keep_evidence is None
        assert parser.parse_args(base+['--keep-evidence']).keep_evidence is True
        assert parser.parse_args(base+['--no-keep-evidence']).keep_evidence is False

def test_reanalysis_current_file_then_set_then_cli(tmp_path,monkeypatch):
    p=parent(tmp_path/'parent',True);seen={};monkeypatch.setattr(reanalysis,'run_audit',lambda **kw:seen.update(kw))
    config=tmp_path/'current.yaml';config.write_text('outputs: {keep_evidence: false}')
    reanalysis.run_reanalysis(results=p,output=tmp_path/'child',config_path=config,sets=['outputs.keep_evidence=true'],keep_evidence=False)
    assert not seen['keep_evidence'] and not yaml.safe_load(Path(seen['config_path']).read_text())['outputs']['keep_evidence']

def test_undeletable_owned_bam_is_not_cleanup_success(tmp_path,monkeypatch):
    out,src,q=finalize_fixture(tmp_path);tmp=out/'evidence/strict_human.selected.unsorted.bam';tmp.write_bytes(b'ARTIFICIAL')
    unlink=Path.unlink
    def denied(self,*a,**kw):
        if self==tmp:raise PermissionError('ARTIFICIAL_DELETE_FAILURE')
        return unlink(self,*a,**kw)
    monkeypatch.setattr(Path,'unlink',denied)
    state=clean(out,False,protected_paths=[src,q])
    assert state['cleanup_status']=='FAILED' and state['cleanup_errors'] and tmp.exists()


def test_real_samtools_sort_failure_removes_actual_bam(tmp_path):
    import shutil
    executable=shutil.which('samtools')
    if not executable:pytest.skip('samtools unavailable')
    out,src,q=finalize_fixture(tmp_path)
    sam=tmp_path/'ARTIFICIAL.sam';sam.write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nARTIFICIAL\t0\tchr1\t1\t60\t100M\t*\t0\t0\t'+'A'*100+'\t'+'I'*100+'\tNM:i:0\n')
    subprocess.run([executable,'view','-b','-o',str(src),str(sam)],check=True,capture_output=True)
    before=src.read_bytes()
    wrapper=tmp_path/'samtools-failing-sort'
    wrapper.write_text('#!'+sys.executable+'\nimport sys,subprocess,pathlib\na=sys.argv[1:]\nif a[0]=="sort":\n pathlib.Path(a[a.index("-o")+1]).write_bytes(b"ARTIFICIAL_PARTIAL_OUTPUT")\n a += ["-m","1K"]\np=subprocess.run(['+repr(executable)+',*a])\nraise SystemExit(p.returncode)\n')
    wrapper.chmod(0o700);dest=out/'evidence/strict_human.bam'
    with pytest.raises(subprocess.CalledProcessError):
        finalize_strict_human_bam(captured_bam=src,qnames_file=q,output_bam=dest,samtools=str(wrapper))
    result=clean(out,False,protected_paths=[src,q])
    assert result['cleanup_status']=='COMPLETE' and not list((out/'evidence').glob('*.bam')) if (out/'evidence').exists() else result['cleanup_status']=='COMPLETE'
    assert src.read_bytes()==before
    assert not dest.with_suffix('.selected.unsorted.bam').exists()

def test_ownership_registered_before_sensitive_creation(tmp_path):
    out,src,q=finalize_fixture(tmp_path);dest=out/'evidence/strict_human.bam'
    def stop(args,**kw):
        record=json.loads((out/'owned_sensitive_intermediates.json').read_text())
        assert {'strict_selection','sort_workspace'} <= {e['kind'] for e in record['entries']}
        assert not dest.with_suffix('.selected.unsorted.bam').exists()
        raise OSError('ARTIFICIAL_STOP')
    with patch('srafeguard.alignment_evidence.subprocess.run',side_effect=stop):
        with pytest.raises(OSError):finalize_strict_human_bam(captured_bam=src,qnames_file=q,output_bam=dest,samtools='ARTIFICIAL')

def test_registered_scratch_is_removed_even_when_final_evidence_retained(tmp_path):
    from srafeguard.owned_temporary import register
    import tempfile
    out=tmp_path/'out';out.mkdir();root=Path(tempfile.mkdtemp(prefix='hla_input_',dir=tmp_path));(root/'ARTIFICIAL.fastq').write_text('ARTIFICIAL')
    neighbor=tmp_path/'user';neighbor.mkdir();(neighbor/'note').write_text('USER')
    register(out,root,'input_workspace');result=clean(out,True)
    assert result['cleanup_status']=='COMPLETE' and not root.exists() and (neighbor/'note').exists()

def test_registered_scratch_deletion_failure_reported(tmp_path,monkeypatch):
    from srafeguard.owned_temporary import register
    import tempfile,shutil
    out=tmp_path/'out';out.mkdir();root=Path(tempfile.mkdtemp(prefix='hla_',dir=tmp_path));register(out,root,'mapping_workspace')
    original=shutil.rmtree
    def fail(path,*args,**kw):
        if Path(path)==root:raise PermissionError('PRIVATE_DETAIL_MUST_NOT_APPEAR')
        return original(path,*args,**kw)
    monkeypatch.setattr(shutil,'rmtree',fail);result=clean(out,False)
    assert result['cleanup_status']=='FAILED' and root.exists()
    assert 'PRIVATE_DETAIL' not in json.dumps(result)

def test_registered_scratch_cannot_remove_original_input(tmp_path):
    from srafeguard.owned_temporary import register
    import tempfile
    out=tmp_path/'out';out.mkdir();root=Path(tempfile.mkdtemp(prefix='hla_input_',dir=tmp_path));original=root/'original.fastq';original.write_text('ORIGINAL')
    register(out,root,'input_workspace');result=clean(out,False,protected_paths=[original])
    assert result['cleanup_status']=='FAILED' and original.read_text()=='ORIGINAL'

def test_registration_cannot_claim_arbitrary_user_directory(tmp_path):
    from srafeguard.owned_temporary import register
    out=tmp_path/'out';out.mkdir();user=tmp_path/'user_data';user.mkdir()
    with pytest.raises(ValueError):register(out,user,'mapping_workspace')

@pytest.mark.parametrize('keep',[False,True])
def test_resolved_cli_policy_written_to_effective_metadata(tmp_path,monkeypatch,keep):
    from test_reporting_compatibility import fake_run
    runner=fake_run(monkeypatch,tmp_path)
    runner.run_audit(sample='ARTIFICIAL',output=tmp_path/'out',keep_evidence=keep)
    out=tmp_path/'out';cfg=yaml.safe_load((out/'effective_config.yaml').read_text());meta=json.loads((out/'run_metadata.json').read_text())
    assert cfg['outputs']['keep_evidence']==keep and meta['effective_parameters']['outputs']['keep_evidence']==keep
    assert meta['evidence_retention_requested']==keep
    assert (out/'logs/human_minimap2.log').exists()==keep
