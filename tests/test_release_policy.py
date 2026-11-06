"""Artificial examples only; no reference-derived human allele profile."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from srafeguard.population import score_population_observations
from srafeguard.runtime_config import load_effective_config
from srafeguard.release_policy import annotate,clean,validate_fastqs,validate_sample,validate_architecture

def observations(n=20):
    return pd.DataFrame([dict(Chromosome='chr1' if i%2 else 'chr2',Position=i*1000000+1,REF='A',ALT='C',Observed_state='ALT',AF_AFR=.5,AF_EUR=.5,AF_SAS=.5,AF_EAS=.5,AF_AMR=.5) for i in range(n)])

def test_support_is_for_full_data_winner_not_max_fraction():
    r=score_population_observations(observations(),bootstraps=100)
    assert max(r['Bootstrap_best_'+p] for p in ['AFR','AMR','EAS','EUR','SAS'])==1
    assert r['Best_bootstrap_support']==r['Bootstrap_best_'+r['Best_population']]
    assert r['Best_bootstrap_support']==0 and not r['Population_affinity_recoverable']

def test_true_zero_missing_and_chromosome_blocks():
    x=observations();x.loc[0,'AF_AFR']=0;x.loc[1,'AF_SAS']=np.nan
    r=score_population_observations(x,bootstraps=100)
    assert r['Sites_used']==19 and not r['Population_affinity_recoverable']
    assert np.isfinite(r['LL_AFR'])
    y=observations(2);y['Position']=1
    assert score_population_observations(y)['Distinct_1Mb_blocks']==2

def test_bootstrap_reproducible_and_eligibility():
    x=observations(25);x['AF_AMR']=.95
    a=score_population_observations(x,seed=31,bootstraps=1000);b=score_population_observations(x,seed=31,bootstraps=1000)
    assert a==b and a['Population_affinity_recoverable']
    assert not score_population_observations(x.iloc[:19],seed=31)['Population_affinity_recoverable']

def report(tmp_path,layout='SINGLE',enabled=True):
    pd.DataFrame([dict(Layout=layout,Population_sites_used=0,Population_affinity_recoverable=False,ClinVar_validated_P_LP_ALT_sites=0,ClinVar_validated_P_LP_supporting_reads=0,ClinVar_ALT_sites_observed=0)]).to_csv(tmp_path/'privacy_summary.tsv',sep='\t',index=False)
    cfg,_=load_effective_config();cfg['clinvar']['enabled']=enabled
    annotate(tmp_path,cfg)
    return pd.read_csv(tmp_path/'privacy_summary.tsv',sep='\t',keep_default_na=False)

def test_single_pair_support_na_is_not_zero(tmp_path):
    x=report(tmp_path);assert x.loc[0,'Pair_supported_P_LP_ALT_sites']=='N/A'
    assert x.loc[0,'Pair_support_applicability']=='N/A_SINGLE_END'
    assert x.loc[0,'Single_read_placement_supported_P_LP_ALT_sites']==0
    assert not any('validated' in c.lower() for c in x)

def test_paired_evaluated_zero_is_zero(tmp_path):
    x=report(tmp_path,'PAIRED');assert x.loc[0,'Pair_supported_P_LP_ALT_sites']==0
    assert x.loc[0,'Pair_support_evaluated']

def test_disabled_clinvar_not_a_negative(tmp_path):
    x=report(tmp_path,'PAIRED',False);assert x.loc[0,'Pair_supported_P_LP_ALT_sites']=='N/A'
    assert x.loc[0,'ClinVar_ALT_sites_observed']=='N/A'

def test_keep_and_clean_retention(tmp_path):
    (tmp_path/'evidence').mkdir();(tmp_path/'evidence/strict_human.bam').write_bytes(b'ARTIFICIAL')
    (tmp_path/'population_marker_observations.tsv.gz').write_bytes(b'ARTIFICIAL')
    clean(tmp_path,True);assert (tmp_path/'evidence/strict_human.bam').exists()
    clean(tmp_path,False);assert not (tmp_path/'evidence').exists()
    assert not (tmp_path/'population_marker_observations.tsv.gz').exists()

def test_installed_defaults_and_config_validation(tmp_path):
    cfg,_=load_effective_config();assert cfg['population']['min_pairwise_af_delta']==.3
    assert cfg['population']['bootstraps']==1000 and not cfg['outputs']['keep_evidence']
    for text in ['strict_human: {min_mapq: 20}','population: {populations: [AFR, EAS]}','population: {af_floor: 0}','population: {bootstrap: {replicates: 100}}']:
        p=tmp_path/'bad.yaml';p.write_text(text)
        with pytest.raises(ValueError):load_effective_config(str(p))

def test_streamed_fastq_qc_rejects_orphans_and_truncation(tmp_path):
    a=tmp_path/'a.fq';b=tmp_path/'b.fq'
    a.write_text('@ARTIFICIAL\nACGT\n+\nIIII\n');b.write_text('@ARTIFICIAL\nACGT\n+\nIIII\n');validate_fastqs(a,b)
    b.write_text('@DIFFERENT/2\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError):validate_fastqs(a,b)
    a.write_text('@ARTIFICIAL\nACGT\n+\n')
    with pytest.raises(ValueError):validate_fastqs(a)

def test_unsafe_labels_and_long_reads_rejected():
    with pytest.raises(ValueError):validate_sample('../outside')
    with pytest.raises(ValueError):validate_architecture('OXFORD_NANOPORE')

def test_nonstandard_paired_sra_indices_fail_before_extraction(tmp_path,monkeypatch):
    from srafeguard import sra_input as si
    from types import SimpleNamespace
    monkeypatch.setattr(si,'_tool',lambda *args:'ARTIFICIAL_TOOL')
    monkeypatch.setattr(si.subprocess,'run',lambda *args,**kwargs:SimpleNamespace(returncode=0))
    monkeypatch.setattr(si,'detect_sra_metadata',lambda *args,**kwargs:dict(layout='PAIRED',biological_read_indices=[1,3]))
    with pytest.raises(ValueError,match='biological read indices'):
        si.extract_and_normalize_sra('ARTIFICIAL.sra','ARTIFICIAL',tmp_path,{})

def test_reanalysis_forwards_recorded_architecture(tmp_path,monkeypatch):
    from srafeguard import reanalysis as rea
    import yaml
    parent=tmp_path/'parent';parent.mkdir()
    (parent/'run_metadata.json').write_text(json.dumps(dict(sample='ARTIFICIAL',inputs=dict(fastq1='ARTIFICIAL.fastq',layout='single',platform='ION_TORRENT'))))
    (parent/'effective_config.yaml').write_text(yaml.safe_dump(dict(outputs=dict(keep_evidence=False))))
    seen={}
    monkeypatch.setattr(rea,'run_audit',lambda **kwargs:seen.update(kwargs))
    rea.run_reanalysis(results=parent,output=tmp_path/'child')
    assert seen['layout']=='single' and seen['platform']=='ION_TORRENT'

def test_reference_builder_restores_zero_and_excludes_missing(tmp_path,monkeypatch):
    from srafeguard import reference_prepare as rp
    import shutil,subprocess,gzip
    if not shutil.which('bcftools'):pytest.skip('bcftools unavailable')
    fields='\n'.join(f'##INFO=<ID={kind}_{p}_unrel,Number=1,Type={"Float" if kind=="AF" else "Integer"},Description="Artificial reference test">' for kind in ['AF','AC','AN'] for p in rp.COLPOPS)
    def info(missing=False):
        parts=[]
        for p in rp.COLPOPS:
            af,ac,an=('.','0','10') if p=='AFR' else ('.','0','0') if missing and p=='AMR' else ('.9','9','10') if p=='AMR' else ('.1','1','10')
            parts.extend([f'AF_{p}_unrel={af}',f'AC_{p}_unrel={ac}',f'AN_{p}_unrel={an}'])
        return ';'.join(parts)
    src=tmp_path/'1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.withafinfo.vcf.gz'
    with gzip.open(src,'wt') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'+fields+'\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        f.write('chr1\t10\tARTIFICIAL\tA\tC\t.\tPASS\t'+info()+'\n')
        f.write('chr1\t20\tARTIFICIAL_MISSING\tA\tC\t.\tPASS\t'+info(True)+'\n')
    monkeypatch.setattr(rp,'DEST',tmp_path);monkeypatch.setattr(rp,'OUT',tmp_path/'out');monkeypatch.setattr(rp,'BCF',Path(shutil.which('bcftools')))
    r=rp.reconstruct(1,'five');assert r['markers']==1 and r['missing_frequency_exclusions']==1
    x=pd.read_csv(tmp_path/'out/five/chr1.tsv.gz',sep='\t');assert x.AF_AFR.tolist()==[0]

def test_provenance_preserves_categorical_architecture(tmp_path,monkeypatch):
    from srafeguard import provenance as pr
    monkeypatch.setattr(pr,'tool_versions',lambda cfg:{})
    monkeypatch.setattr(pr,'resource_hashes',lambda *args:{})
    p=pr.start_run_record(output=tmp_path,sample='ARTIFICIAL',cfg={},config_path=None,inputs=dict(fastq1='artificial.fastq',layout='auto',platform='ION_TORRENT'))
    data=json.loads(p.read_text())
    assert data['inputs']['layout']=='auto'
    assert data['inputs']['platform']=='ION_TORRENT'
    assert Path(data['inputs']['fastq1']).is_absolute()
