"""Portable reference preparation; five-population extraction is the frozen producer.

No reference data are bundled or downloaded by this module. Original full-source
controls matched the frozen four-group panel before the final five-group build.
"""
import argparse,csv,gzip,hashlib,itertools,json,math,pathlib,shutil,struct,subprocess
POPS4=['AFR','EUR','SAS','EAS'];POPS5=['AFR','AMR','EAS','EUR','SAS'];COLPOPS=POPS4+['AMR']
FIELDS=['Chromosome','Position','REF','ALT']+[f'AF_{p}' for p in COLPOPS]+['Max_pairwise_AF_delta','Most_separated_populations']
DEST=OUT=BCF=None

def digest(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def output_gzip(path):
    raw=path.open('wb'); gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)
    import io
    return io.TextIOWrapper(gz,encoding='utf-8',newline='')

def delta(freq,pops):
    pairs=[(abs(freq[a]-freq[b]),f'{a}-{b}') for a,b in itertools.combinations(pops,2)]
    return max(pairs,key=lambda p:p[0])

def number(x):
    if x in {'','.','NA','nan','None'}: return None
    v=float(x)
    if not math.isfinite(v) or not 0<=v<=1: raise ValueError('Invalid AF')
    return v

def reconstruct(c,mode):
    src=DEST/f'1KGP.CHM13v2.0.chr{c}.recalibrated.snp_indel.pass.withafinfo.vcf.gz'
    if not src.exists(): raise FileNotFoundError(src)
    out=OUT/mode/f'chr{c}.tsv.gz'; out.parent.mkdir(parents=True,exist_ok=True)
    fields=['%CHROM','%POS','%REF','%ALT']+[f'%INFO/AF_{p}_unrel' for p in COLPOPS]
    fields += [f'%INFO/{kind}_{p}_unrel' for kind in ['AC','AN'] for p in COLPOPS]
    cmd=[str(BCF),'query','-f','\t'.join(fields)+'\n',str(src)]
    stat={'chromosome':f'chr{c}','source_sha256':digest(src),'mode':mode,'command':cmd,'records':0,'biallelic_acgt_snv':0,'eligible_complete_frequency_universe':0,'missing_frequency_exclusions':0,'recovered_zero_frequencies':{p:0 for p in COLPOPS},'qualified_before_ambiguity':0,'ambiguous_coordinates':0,'markers':0,'threshold_counts':{str(t):0 for t in [.1,.2,.3]}}
    stderr_path=out.with_suffix('.stderr.txt')
    with stderr_path.open('w') as err, output_gzip(out) as dest:
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1024*1024)
        writer=csv.writer(dest,delimiter='\t',lineterminator='\n'); writer.writerow(FIELDS)
        pending=[]; previous=None; lastpos=0
        def flush():
            if len(pending)>1: stat['ambiguous_coordinates']+=1
            elif pending:
                row=pending[0]; writer.writerow(row); stat['markers']+=1
                for t in [.1,.2,.3]: stat['threshold_counts'][str(t)]+=row[-2]>=t
        try:
            for line in proc.stdout:
                f=line.rstrip('\n').split('\t')
                if len(f)!=19: raise ValueError(f'query columns {len(f)} !=19')
                chrom,pos,ref,alt=f[:4]; pos=int(pos); stat['records']+=1
                if chrom!=f'chr{c}' or pos<lastpos: raise ValueError('Unexpected chromosome or unsorted coordinate')
                lastpos=pos
                if previous!=(chrom,pos): flush(); pending=[]; previous=(chrom,pos)
                ref=ref.upper(); alt=alt.upper()
                if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT': continue
                stat['biallelic_acgt_snv']+=1
                af={p:number(f[4+i]) for i,p in enumerate(COLPOPS)}
                if mode=='five':
                    for i,p in enumerate(COLPOPS):
                        ac,an=f[9+i],f[14+i]
                        if af[p] is None and ac=='0' and an not in {'.',''} and int(an)>0:
                            af[p]=0.; stat['recovered_zero_frequencies'][p]+=1
                        if af[p] is not None:
                            if ac in {'.',''} or an in {'.',''} or int(an)<=0 or not 0<=int(ac)<=int(an): raise ValueError('Invalid authoritative AC/AN')
                            # The source Perl formatter emits 8 decimal places;
                            # VCF Float is float32 and bcftools query prints 6
                            # significant digits. Check that exact documented
                            # conversion, rather than a guessed error tolerance.
                            rounded8=float(f'{int(ac)/int(an):.8f}')
                            single=struct.unpack('f',struct.pack('f',rounded8))[0]
                            expected=float(f'{single:.6g}')
                            if af[p]!=expected:raise ValueError('AF does not match authoritative counts and exact source/query formatting')
                pops=POPS4 if mode=='control' else POPS5
                if any(af[p] is None for p in pops): stat['missing_frequency_exclusions']+=1; continue
                stat['eligible_complete_frequency_universe']+=1
                d,pair=delta(af,pops)
                if d<.10: continue
                stat['qualified_before_ambiguity']+=1
                pending.append([chrom,pos,ref,alt]+[af[p] if af[p] is not None else '' for p in COLPOPS]+[d,pair])
            flush()
        except BaseException:
            proc.terminate(); proc.wait(); raise
        if proc.wait()!=0: raise RuntimeError(f'bcftools failed; see {stderr_path}')
    stat['output_sha256']=digest(out)
    out.with_suffix('.receipt.json').write_text(json.dumps(stat,indent=2)+'\n')
    print(json.dumps({k:stat[k] for k in ['chromosome','mode','records','markers','ambiguous_coordinates']}),flush=True)
    return stat

def main():
    global DEST,OUT,BCF
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',required=True,type=pathlib.Path);p.add_argument('--output-dir',required=True,type=pathlib.Path);a=p.parse_args()
    DEST=a.input_dir.resolve();OUT=a.output_dir.resolve();BCF=pathlib.Path(shutil.which('bcftools') or 'bcftools')
    if OUT.exists():raise FileExistsError('Use a new output directory')
    manifest=json.loads(pathlib.Path(__file__).with_name('data').joinpath('reference_manifest.json').read_text())
    for item in manifest['chromosomes']:
        path=DEST/item['file']
        if digest(path)!=item['sha256']:raise ValueError('Upstream chromosome checksum differs: '+item['file'])
    OUT.mkdir(parents=True)
    chromosomes=sorted(range(1,23),key=lambda c:f'chr{c}')
    stats=[reconstruct(c,'five') for c in chromosomes]
    target=OUT/'five/chm13_1kgp_five_population_delta010.tsv.gz'
    with output_gzip(target) as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(FIELDS)
        for c in chromosomes:
            with gzip.open(OUT/'five'/f'chr{c}.tsv.gz','rt') as source:
                reader=csv.reader(source,delimiter='\t');assert next(reader)==FIELDS
                writer.writerows(reader)
    actual=digest(target)
    receipt=dict(status='PASS' if actual==manifest['resource_sha256'] else 'FAIL_STOP',sha256=actual,chromosomes=stats)
    (OUT/'reference_build_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if receipt['status']!='PASS':raise RuntimeError('Output does not match frozen final resource; do not use it')
    print(json.dumps(dict(status='PASS',resource=str(target),sha256=actual)))

if __name__=='__main__':main()
