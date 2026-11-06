"""Installed-resource configuration; no research-tree auto-discovery."""
from pathlib import Path
import os,tempfile,copy
import yaml
DATA=Path(__file__).with_name('data')

def merge(base,over):
    for k,v in over.items():
        if isinstance(v,dict) and isinstance(base.get(k),dict):merge(base[k],v)
        else:base[k]=copy.deepcopy(v)
    return base

def load(explicit=None):
    cfg=yaml.safe_load((DATA/'default.yaml').read_text())
    path=Path(explicit).expanduser().resolve() if explicit else None
    if path:
        if not path.is_file():raise FileNotFoundError(path)
        over=yaml.safe_load(path.read_text()) or {}
        if not isinstance(over,dict):raise ValueError('configuration must be a mapping')
        merge(cfg,over)
    base=path.parent if path else Path.cwd()
    for k,v in cfg.get('resources',{}).items():
        if v=='@bundled_mask':cfg['resources'][k]=str(DATA/'extreme_hotspot_exclusion_mask.bed')
        elif isinstance(v,str) and v:
            p=Path(v).expanduser();cfg['resources'][k]=str((base/p).resolve() if not p.is_absolute() else p)
    rt=cfg.setdefault('runtime',{})
    if rt.get('keep_intermediates',False):raise ValueError('use the documented keep-evidence control; retained scratch intermediates are not supported by the public cleanup policy')
    work=rt.get('temp_dir') or rt.get('tempdir') or str(Path(tempfile.gettempdir())/('srafeguard-'+str(os.getuid())))
    work=Path(work).expanduser();work=(base/work).resolve() if not work.is_absolute() else work
    rt['temp_dir']=rt['tempdir']=str(work)
    fixed={'min_mapq':30,'min_identity':0.98,'min_query_coverage':0.8,'require_proper_pair':True,'require_both_mates':True,'primary_only':True,'deduplicate_physical_fragments':True,'mapper_preset':'sr'}
    if any(cfg.get('strict_human',{}).get(k)!=v for k,v in fixed.items()):
        raise ValueError('strict caller thresholds are frozen; these configuration keys cannot change the scientific implementation')
    pop=cfg.get('population',{})
    if cfg.get('release',{}).get('reference_profile')=='FINAL_FIVE_POPULATION' and float(pop.get('min_pairwise_af_delta',.3))<.1:
        raise ValueError('the final resource starts at Delta AF .10; lower thresholds require a separately validated complete CUSTOM resource')
    for unused in ['superpopulations','bootstrap','remap','interpretation']:
        if unused in pop:raise ValueError('use documented flat population configuration fields, not '+unused)
    if len(pop.get('populations',[]))!=5 or set(pop.get('populations',[]))!={'AFR','AMR','EAS','EUR','SAS'}:raise ValueError('five reference groups required')
    if not 0<float(pop.get('af_floor',.005))<.5:raise ValueError('invalid AF floor')
    if not 0<=float(pop.get('min_pairwise_af_delta',.3))<=1:raise ValueError('invalid Delta AF threshold')
    if int(pop.get('bootstraps',0))<0 or int(pop.get('informative_min_sites',0))<1 or int(pop.get('block_size_bp',0))<1:raise ValueError('invalid population counts/block size')
    if not 0<=float(pop.get('informative_min_bootstrap',.9))<=1:raise ValueError('invalid support threshold')
    if cfg.get('input',{}).get('assay','genomic')!='genomic':raise ValueError('splice-aware transcriptomic assay is not implemented')
    reporting=cfg.get('reporting',{}) or {}
    if reporting.get('schema','public-2') not in ('public-2','legacy-1'):raise ValueError('reporting.schema must be public-2 or legacy-1')
    if not isinstance(reporting.get('legacy_export',False),bool):raise ValueError('reporting.legacy_export must be boolean')
    return cfg
