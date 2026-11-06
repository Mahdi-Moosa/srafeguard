from srafeguard.batch import load_manifest
from srafeguard.provenance import _sha256
from srafeguard.reanalysis import _apply_common_overrides
from srafeguard.runner import _infer_sample

def test_manifest_loader(tmp_path):
    p=tmp_path/"manifest.tsv"
    p.write_text("sample\tsra\tfastq1\tfastq2\treject_source_sra\nA\t/a.sra\t\t\t\nB\t\t/b1.fq.gz\t/b2.fq.gz\t/b.sra\n")
    x=load_manifest(p)
    assert list(x["sample"])==["A","B"]

def test_reanalysis_changes_only_explicit_fields():
    cfg={"population":{"min_mapq":30,"min_bq":30,"other":7},
         "clinvar":{"min_mapq":30,"min_bq":30}}
    out,overrides=_apply_common_overrides(cfg,min_bq=20,sets=["clinvar.chromosome_scope=autosomes"])
    assert out["population"]["min_mapq"]==30
    assert out["population"]["min_bq"]==20
    assert out["population"]["other"]==7
    assert out["clinvar"]["min_bq"]==20
    assert out["clinvar"]["chromosome_scope"]=="autosomes"
    assert "min_mapq" not in overrides

def test_sha256(tmp_path):
    p=tmp_path/"x"; p.write_bytes(b"abc")
    assert _sha256(p)=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

def test_sample_inference():
    assert _infer_sample(fastq1="ABC_1.fastq.gz")=="ABC"
