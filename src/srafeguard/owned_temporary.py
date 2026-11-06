"""Ownership of temporary evidence, registered before sensitive data are written."""
from pathlib import Path
import json,os,re,shutil
REGISTRY='owned_sensitive_intermediates.json'
PATTERNS={
    'strict_selection':r'strict_human\.selected\.unsorted\.bam',
    'captured_mapping':r'[A-Za-z0-9_.-]+\.original_human_mapping\.bam',
    'input_workspace':r'hla_input_[A-Za-z0-9_.-]+',
    'mapping_workspace':r'hla_[A-Za-z0-9_.-]+',
    'sort_workspace':r'srafeguard_sort_[A-Za-z0-9_.-]+',
    'clinical_workspace':r'srafeguard_clinvar_[A-Za-z0-9_.-]+',
}
DIRECTORIES={'input_workspace','mapping_workspace','sort_workspace','clinical_workspace'}

def _validate(path,kind):
    if kind not in PATTERNS or not re.fullmatch(PATTERNS[kind],path.name):
        raise ValueError('Unrecognized temporary ownership record')
    if path==Path(path.anchor) or path==Path.home():raise ValueError('Unsafe temporary ownership record')

def _load(output):
    p=Path(output)/REGISTRY
    if p.is_symlink():raise ValueError('Ownership registry must not be a symlink')
    if not p.exists():return dict(schema='srafeguard-owned-temporary-1',output=str(Path(output).resolve()),entries=[])
    d=json.loads(p.read_text())
    if d.get('schema')!='srafeguard-owned-temporary-1' or d.get('output')!=str(Path(output).resolve()):raise ValueError('Invalid temporary ownership registry')
    return d

def register(output,path,kind):
    """Only producer-controlled names and freshly created scratch directories qualify."""
    output=Path(output);path=Path(path).absolute();_validate(path,kind)
    if path.is_symlink():raise ValueError('Temporary evidence must not be a symlink')
    d=_load(output)
    if not any(e['path']==str(path) for e in d['entries']):
        item=dict(path=str(path),kind=kind)
        if kind in DIRECTORIES:
            st=path.stat();item.update(device=st.st_dev,inode=st.st_ino)
        d['entries'].append(item)
    p=output/REGISTRY;tmp=output/(REGISTRY+'.tmp')
    with tmp.open('w') as f:json.dump(d,f,indent=2)
    os.replace(tmp,p)

def clean_registered(output,protected_paths=()):
    """Clean exclusively registered ephemeral work, independently of retained final evidence."""
    errors=[];protected=[Path(p).resolve() for p in protected_paths if p]
    try:d=_load(output)
    except Exception as exc:return [dict(file=REGISTRY,error_type=type(exc).__name__)]
    for e in d['entries']:
        try:
            p=Path(e['path']);kind=e['kind'];_validate(p,kind)
            try:st=p.lstat()
            except FileNotFoundError:continue
            if any(p.resolve()==q or (kind in DIRECTORIES and q.is_relative_to(p.resolve())) for q in protected):raise ValueError('Original input is protected')
            if p.is_symlink():p.unlink()
            elif kind in DIRECTORIES:
                if (st.st_dev,st.st_ino)!=(e.get('device'),e.get('inode')):raise ValueError('Temporary directory identity changed')
                shutil.rmtree(p)
            else:p.unlink()
            try:p.lstat()
            except FileNotFoundError:continue
            raise OSError('Owned temporary remains')
        except Exception as exc:errors.append(dict(file='TEMPORARY:'+e.get('kind','UNKNOWN'),error_type=type(exc).__name__))
    return errors
