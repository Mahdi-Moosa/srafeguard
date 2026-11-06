"""Generate and audit only artificial references, allele frequencies and VCF annotations."""
from pathlib import Path
import json,random,csv,subprocess,os,time,shutil
import yaml,pandas as pd
import argparse,sys
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workdir',required=True);args=parser.parse_args()
D=Path(args.workdir).resolve()
if D.exists():raise FileExistsError('Choose a new work directory outside the source tree')
D.mkdir(mode=0o700,parents=True);B=D;E=Path(sys.executable).parent.parent;(B/'review').mkdir(mode=0o700);os.umask(0o077)
env=os.environ.copy();env.pop('PYTHONPATH',None);env['PATH']=str(E/'bin')+os.pathsep+os.environ.get('PATH','');env['PYTHONDONTWRITEBYTECODE']='1';env['XDG_CACHE_HOME']=str(B/'cache')
rng=random.Random('ARTIFICIAL_REFERENCE_ONLY_NOT_HUMAN');ref={c:''.join(rng.choice('ACGT') for _ in range(1200000)) for c in ['chr1','chr2']}
fa=D/'ARTIFICIAL_reference.fa'
with fa.open('w') as f:
 for c,s in ref.items():f.write('>'+c+'\n'+s+'\n')
r1=[];r2=[];markers=[];variants=[]
for chrom in ref:
 for block in [0,1]:
  for i in range(8):
   left=block*1000000+1000+i*1000;seq=ref[chrom][left:left+150];pos=left+76;base=ref[chrom][pos-1];alt=next(b for b in 'ACGT' if b!=base);seq=seq[:75]+alt+seq[76:]
   mate=ref[chrom][left+250:left+400].translate(str.maketrans('ACGT','TGCA'))[::-1];name='ARTIFICIAL_PAIR_'+str(len(r1));r1.append((name,seq));r2.append((name,mate))
   markers.append(dict(Chromosome=chrom,Position=pos,REF=base,ALT=alt,Max_pairwise_AF_delta=.95,AF_AFR=0,AF_EUR=.2,AF_SAS=.3,AF_EAS=.05,AF_AMR=.95))
   if len(r1)==1:variants.append((chrom,pos,'ARTIFICIAL_PL',base,alt,'Pathogenic'))
   if len(r1)==2:variants.append((chrom,pos,'ARTIFICIAL_BENIGN',base,alt,'Benign'))
   if len(r1)==3:
    p=left+81;b=ref[chrom][p-1];a=next(x for x in 'ACGT' if x!=b);variants.append((chrom,p,'ARTIFICIAL_REF',b,a,'Likely_pathogenic'))
for name,rows in [('R1',r1),('R2',r2)]:
 with (D/(name+'.fastq')).open('w') as f:
  for q,s in rows:f.write('@'+q+'\n'+s+'\n+\n'+'I'*len(s)+'\n')
pd.DataFrame(markers).to_csv(D/'ARTIFICIAL_markers.tsv.gz',sep='\t',index=False,compression='gzip')
vcf=D/'ARTIFICIAL_annotations.vcf'
with vcf.open('w') as f:
 f.write('##fileformat=VCFv4.2\n##reference=WHOLLY_ARTIFICIAL_NOT_HUMAN\n')
 for c in ref:f.write('##contig=<ID='+c+',length=1200000>\n')
 f.write('##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Artificial test label, not a clinical record">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
 for c,p,i,r,a,label in variants:f.write(f'{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\tCLNSIG={label}\n')
subprocess.run([(shutil.which('minimap2',path=env['PATH']) or 'minimap2'),'-x','sr','-d',str(D/'ARTIFICIAL.sr.mmi'),str(fa)],env=env,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
subprocess.run([(shutil.which('bgzip',path=env['PATH']) or 'bgzip'),'-f',str(vcf)],env=env,check=True);subprocess.run([(shutil.which('tabix',path=env['PATH']) or 'tabix'),'-p','vcf',str(vcf)+'.gz'],env=env,check=True)
cfg=dict(profile=dict(name='CUSTOM_ARTIFICIAL_ONLY'),release=dict(reference_profile='CUSTOM_SYNTHETIC'),resources=dict(human_reference=str(fa),human_minimap2_index=str(D/'ARTIFICIAL.sr.mmi'),population_markers=str(D/'ARTIFICIAL_markers.tsv.gz'),clinvar_vcf=str(vcf)+'.gz',clinvar_vcf_index=str(vcf)+'.gz.tbi'),runtime=dict(temp_dir=str(D/'tmp'),mapper_threads=2,sort_threads=1),outputs=dict(keep_evidence=False))
(D/'synthetic.yaml').write_text(yaml.safe_dump(cfg));checks=[]
def execute(step,args):
 start=time.time();r=subprocess.run([str(E/'bin/srafeguard'),*args],env=env,cwd=D,capture_output=True,text=True,timeout=240)
 (D/(step+'.log')).write_text(r.stdout+r.stderr);checks.append(dict(step=step,arguments=args,returncode=r.returncode,seconds=time.time()-start))
 if r.returncode:raise RuntimeError(step+' failed')
base=['--fastq1',str(D/'R1.fastq'),'--fastq2',str(D/'R2.fastq'),'--sample','ARTIFICIAL_PAIRED','--config',str(D/'synthetic.yaml')]
execute('paired',['audit',*base,'--output',str(D/'paired')])
x=pd.read_csv(D/'paired/privacy_summary.tsv',sep='\t',keep_default_na=False).iloc[0]
assert x['Strict_human_fragments']==32 and x['Human_chromosomes_represented']==2 and x['Population_sites_used']==32 and x['Population_distinct_1Mb_blocks']==4
assert x['Population_best_affinity']=='AMR' and x['Population_best_bootstrap']==1 and x['Population_affinity_recoverable']
assert x['ClinVar_sites_with_HQ_base_observation']==3 and x['ClinVar_ALT_sites_observed']==2 and x['ClinVar_candidate_P_LP_ALT_sites']==1 and x['Pair_supported_P_LP_ALT_sites']==1
assert not (D/'paired/evidence/strict_unique_fragments.tsv').exists() and not (D/'paired/population_marker_observations.tsv.gz').exists()
execute('retained',['audit',*base,'--keep-evidence','--output',str(D/'retained')]);assert (D/'retained/evidence/strict_human.bam').is_file();assert (D/'retained/population_marker_observations.tsv.gz').exists()
y=pd.read_csv(D/'retained/privacy_summary.tsv',sep='\t',keep_default_na=False).iloc[0]
assert x.to_dict()==y.to_dict()
execute('single',['audit','--fastq1',str(D/'R1.fastq'),'--layout','single','--platform','ION_TORRENT','--sample','ARTIFICIAL_SINGLE','--config',str(D/'synthetic.yaml'),'--output',str(D/'single')])
z=pd.read_csv(D/'single/privacy_summary.tsv',sep='\t',keep_default_na=False).iloc[0];assert z['Pair_supported_P_LP_ALT_sites']=='N/A';assert z['Single_read_placement_supported_P_LP_ALT_sites']==1
execute('reanalyze',['reanalyze','--results',str(D/'paired'),'--output',str(D/'reanalyzed')]);a=pd.read_csv(D/'reanalyzed/privacy_summary.tsv',sep='\t',keep_default_na=False).iloc[0];assert a.to_dict()==x.to_dict()
with (D/'batch.tsv').open('w') as f:
 f.write('sample\tfastq1\tfastq2\tlayout\nARTIFICIAL_BATCH\t'+str(D/'R1.fastq')+'\t'+str(D/'R2.fastq')+'\tpaired\n')
execute('batch',['batch','--manifest',str(D/'batch.tsv'),'--output-root',str(D/'batch_results'),'--config',str(D/'synthetic.yaml')]);execute('batch_resume',['batch','--manifest',str(D/'batch.tsv'),'--output-root',str(D/'batch_results'),'--config',str(D/'synthetic.yaml')])
batch=pd.read_csv(D/'batch_results/batch_status.tsv',sep='\t');assert batch.status.tolist()==['skipped_completed']
receipt=dict(status='PASS',data_origin='Wholly artificial generated reference/reads/AFs/VCF labels; no human-derived sequences, alleles or clinical examples.',checks=checks,synthetic_expected=dict(strict_pairs=32,chromosomes=2,clean_markers=32,blocks=4,HQ_annotated_sites=3,ALT_sites=2,candidate_PL_sites=1,pair_supported_sites=1),retention_preserves_numerical_results=True,single_pair_support_NA=True,single_placement_is_separate=True,reanalysis_matches=True,batch_resume_skips=True,default_details_removed=True,historical_clinical_or_affinity_validation_claim=False)
(B/'review/synthetic_integration.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='checks'}))
