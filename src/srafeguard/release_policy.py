"""Public report applicability and explicit evidence-retention policy."""
from pathlib import Path
import json,os,re,hashlib,shutil
import pandas as pd

def validate_sample(sample):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',str(sample)) or str(sample) in {'.','..'}:
        raise ValueError('sample must be a simple filename-safe label')

def validate_architecture(platform):
    if any(x in str(platform).upper() for x in ['NANOPORE','PACBIO','OXFORD','SEQUEL']):
        raise ValueError('long-read architecture is not validated by the short-read callers')

def trusted_resources(cfg):
    if cfg.get('release',{}).get('reference_profile')!='FINAL_FIVE_POPULATION':return
    resources=cfg.get('resources',{})
    checks=[]
    if cfg.get('population',{}).get('enabled',True):checks.append(('population_markers','6e5655a095f1b071915e5c5f527a4a698e96dcf08a6cd58cba141f7a7ecdc112'))
    checks.append(('human_reference','a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3'))
    if cfg.get('clinvar',{}).get('enabled',True):checks.append(('clinvar_vcf','a551d3339915c038aaf261200bbe1dc7d25203a44ea6141eac8efc07f9b66d8c'))
    for key,expected in checks:
        value=resources.get(key)
        if not value or not Path(value).is_file():raise FileNotFoundError('configure verified external resource '+key)
        h=hashlib.sha256()
        with Path(value).open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
        if h.hexdigest()!=expected:raise ValueError('resource differs from finalized five-population profile: '+key)

# Names are relative to the audit output. Unknown files and directories are untouched.
SENSITIVE_FILES = (
    'evidence/strict_human.bam', 'evidence/strict_human.bam.bai',
    'evidence/strict_human.selected.unsorted.bam',
    'evidence/strict_unique_R1.fastq.gz', 'evidence/strict_unique_R2.fastq.gz',
    'evidence/strict_unique_R1.fastq', 'evidence/strict_unique_R2.fastq',
    'evidence/strict_unique.fastq.gz', 'evidence/strict_unique.fastq',
    'evidence/strict_unique_qnames.txt', 'evidence/strict_unique_fragments.tsv',
    'evidence/strict_unique_reads.tsv', 'evidence/exact_bacterial.tmp.bam',
    'evidence/exact_bacterial.tmp.bam.bai', 'evidence/bacterial_minimap2.log',
    'logs/human_minimap2.log', 'population_marker_observations.tsv.gz',
    'clinvar_site_observations_raw.tsv.gz', 'clinvar_site_observations.tsv.gz',
    'clinvar_plp_read_validation.tsv.gz',
)
SCHEMA = 'srafeguard-public-2'
LEGACY_FIELDS = (
    'ClinVar_validated_P_LP_ALT_sites', 'ClinVar_validated_P_LP_supporting_reads',
    'ClinVar_validated_P_LP_supporting_mates', 'P_LP_validated_exposure_detected',
)

def annotate(output,cfg):
    """Wrap the untouched scientific report; legacy exports preserve original values."""
    output=Path(output);p=output/'privacy_summary.tsv'
    if not p.is_file():return
    x=pd.read_csv(p,sep='\t',keep_default_na=False)
    legacy=x.copy(deep=True)
    reporting=cfg.get('reporting',{}) or {}
    if reporting.get('legacy_export',False) or reporting.get('schema')=='legacy-1':
        legacy.to_csv(output/'privacy_summary_legacy.tsv',sep='\t',index=False)
    layout=x['Layout'].astype(str).str.upper()
    x['Strict_units']=x.get('Strict_human_fragments','')
    x['Counting_unit']=x.get('Evidence_unit','')
    enabled=bool(cfg.get('clinvar',{}).get('enabled',True))
    validation=enabled and bool(cfg.get('clinvar',{}).get('validation',{}).get('enabled',True))
    pop=bool(cfg.get('population',{}).get('enabled',True))
    values=x.get('ClinVar_validated_P_LP_ALT_sites',pd.Series(['']*len(x)))
    mates=x.get('ClinVar_validated_P_LP_supporting_mates',pd.Series(['']*len(x)))
    reads=x.get('ClinVar_validated_P_LP_supporting_reads',pd.Series(['']*len(x)))
    x['Pair_support_applicability']=['APPLICABLE' if z=='PAIRED' else 'N/A_SINGLE_END' for z in layout]
    x['Pair_support_evaluated']=validation & layout.eq('PAIRED')
    x['Single_read_placement_evaluated']=validation & layout.eq('SINGLE')
    x['Pair_supported_P_LP_ALT_sites']=[v if z=='PAIRED' and validation else 'N/A' for v,z in zip(values,layout)]
    x['Pair_supported_P_LP_supporting_mates']=[v if z=='PAIRED' and validation else 'N/A' for v,z in zip(mates,layout)]
    x['Single_read_placement_supported_P_LP_ALT_sites']=[v if z=='SINGLE' and validation else 'N/A' for v,z in zip(values,layout)]
    x['Single_read_placement_supporting_reads']=[v if z=='SINGLE' and validation else 'N/A' for v,z in zip(reads,layout)]
    x['ClinVar_evaluation_status']='EVALUATED' if enabled else 'NOT_EVALUATED_DISABLED'
    x['Population_evaluation_status']='EVALUATED' if pop else 'NOT_EVALUATED_DISABLED'
    x['Population_marker_positive']=(pd.to_numeric(x['Population_sites_used'],errors='coerce')>0) if pop else 'N/A'
    x['Population_marker_eligible']=(pd.to_numeric(x['Population_sites_used'],errors='coerce')>=cfg['population']['informative_min_sites']) if pop else 'N/A'
    x=x.drop(columns=[c for c in LEGACY_FIELDS if c in x])
    config_fields={'ClinVar_min_MAPQ','ClinVar_min_BQ','ClinVar_chromosome_scope','ClinVar_apply_hotspot_exclusion_mask','ClinVar_validation_mode','Population_min_MAPQ','Population_min_BQ'}
    if not enabled:
        for col in x:
            if col.startswith(('ClinVar_','Raw_ClinVar_','Hotspot_','P_LP_')) and col not in config_fields|{'ClinVar_evaluation_status'}:x[col]='N/A'
    if not pop:
        for col in x:
            if col.startswith('Population_') and col not in config_fields|{'Population_evaluation_status'}:x[col]='N/A'
    x['Report_schema']=SCHEMA
    if reporting.get('schema','public-2')=='legacy-1':
        x.to_csv(output/'privacy_summary_public.tsv',sep='\t',index=False)
        legacy.to_csv(p,sep='\t',index=False)
    else:x.to_csv(p,sep='\t',index=False)
    if len(layout) and layout.iloc[0]=='SINGLE':
        meta=output/'clinvar_analysis_metadata.json'
        if meta.is_file():
            d=json.loads(meta.read_text());d['validation_mode']='single_read_unambiguous_remap';d.setdefault('validation',{})['rule']='Independent read ALT/locus placement support; physical read-pair confirmation is N/A for SINGLE input.';meta.write_text(json.dumps(d,indent=2)+'\n')

def clean(output,keep,protected_paths=()):
    """Delete only enumerated owned detail files; errors are explicit, never suppressed."""
    from .owned_temporary import clean_registered
    output=Path(output);errors=clean_registered(output,protected_paths);removed=[]
    protected={Path(p).resolve() for p in protected_paths if p}
    if not keep:
        for name in SENSITIVE_FILES:
            p=output/name
            try:
                if any((output/q).is_symlink() for q in p.relative_to(output).parents if str(q)!='.'):
                    raise ValueError('symlink ancestor')
                if p.resolve() in protected:raise ValueError('protected original input')
                try:p.lstat()
                except FileNotFoundError:continue
                p.unlink();removed.append(name)
                try:p.lstat()
                except FileNotFoundError:continue
                raise OSError('Owned sensitive file remains')
            except Exception as exc:errors.append({'file':name,'error_type':type(exc).__name__})
        # Remove only empty directories; unrelated user files remain.
        for name in ('evidence','logs'):
            p=output/name
            if p.is_dir() and not p.is_symlink():
                try:
                    if not any(p.iterdir()):p.rmdir()
                except OSError as exc:errors.append({'file':name,'error_type':type(exc).__name__})
    result=dict(schema=SCHEMA,keep_evidence=bool(keep),cleanup_status='FAILED' if errors else 'COMPLETE',
                evidence_retention_status='PARTIAL_CLEANUP_FAILURE' if errors else 'RETAINED_BY_REQUEST' if keep else 'SENSITIVE_DETAILS_REMOVED',
                removed_files=removed,cleanup_errors=errors,summary_is_not_anonymous=True)
    try:(output/'retention_policy.json').write_text(json.dumps(result,indent=2)+'\n')
    except Exception as exc:
        result['cleanup_status']='FAILED';result['cleanup_errors'].append({'file':'retention_policy.json','error_type':type(exc).__name__})
        result['evidence_retention_status']='PARTIAL_CLEANUP_FAILURE'
    return result

def diagnostic(output,phase,error_type):
    """A minimal diagnostic without raw exception text, sequences, QNAMEs or commands."""
    (Path(output)/'diagnostic_summary.json').write_text(json.dumps(dict(phase=phase,error_type=error_type,details='Consult run status and configuration; raw tool logs removed unless evidence retention was requested.'))+'\n')


def validate_fastqs(fastq1,fastq2=None):
    if fastq1 is None:return
    import gzip,itertools
    def records(path):
        opener=gzip.open if str(path).endswith('.gz') else open
        with opener(path,'rt') as f:
            while True:
                h=f.readline()
                if not h:break
                seq,plus,qual=f.readline(),f.readline(),f.readline()
                if not h.startswith('@') or not plus.startswith('+') or not qual or len(seq.rstrip())!=len(qual.rstrip()):raise ValueError('invalid or incomplete FASTQ record')
                name=h[1:].split()[0]
                if fastq2 and name.endswith(('/1','/2')):raise ValueError('Normalize mate suffixes externally; paired strict-read exports require identical identifiers')
                yield name
    if fastq2:
        for a,b in itertools.zip_longest(records(fastq1),records(fastq2)):
            if a is None or b is None or a!=b:raise ValueError('paired FASTQ records must have synchronized matching identifiers and counts')
    else:
        for _ in records(fastq1):pass
