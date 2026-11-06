# Reference preparation and licenses

References are external except the small study-curated BED exclusion mask. No
human FASTA/index, 1KGP dataset, ClinVar VCF or actual genetic evidence is bundled.
Do not infer data redistribution rights from the wrapper MIT license. Review each
provider's resource terms independently. Software dependencies retain their own
licenses; no third-party binary/source is vendored.

## Human analysis set

Recorded source:
`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz`

This includes CHM13v2.0 autosomes/X, the completed Y with masked pseudoautosomal
regions, and rCRS mitochondrial sequence; it does not mask the entire Y chromosome.
Upstream compressed MD5: `b17cd2bdd8eb5544e195d6a8882e417b`.
Uncompressed FASTA SHA-256:
`a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3`.

```bash
# Download/decompress the named provider resource outside this repository.
# Verify compressed MD5, then uncompressed SHA-256 against the values above.
minimap2 -x sr -d /path/to/external/chm13.sr.mmi /path/to/external/chm13.fa
samtools faidx /path/to/external/chm13.fa
```

Save executable version and index SHA-256. The retained historical short-read index
hash is `2158acb5aef0f7b9e8c09fcab1ac7e3c9a7281817f0a9a081380074580a72555`.
A rebuilt index can differ across tool versions; retain its exact lineage to the
verified FASTA. Do not substitute asm/map-ont indexes or another human reference.

## Complete five-population marker resource

The final resource is `chm13_1kgp_five_population_delta010.tsv.gz`, SHA-256:
`6e5655a095f1b071915e5c5f527a4a698e96dcf08a6cd58cba141f7a7ecdc112`.
It has 6,947,747/.10, 4,013,865/.20 and 2,156,198/.30 eligible markers. SRAFEGUARD
uses .30 by default, configurable using the same complete table at .10 or higher.

Source release base:
`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/unrelated_samples_2504/allele_freq/`

Autosomal names:
`1KGP.CHM13v2.0.chrN.recalibrated.snp_indel.pass.withafinfo.vcf.gz`, N=1..22.
The complete source is about8.83GiB compressed. Obtain the sites/allele-frequency
release, not a different coordinate system or reference panel. Provider URLs and
SHA-256 of all22 files are installed in `data/reference_manifest.json`.

After retrieval outside this repository, the included portable preparation entry
point verifies all22 upstream checksums, runs the frozen five-population extraction
function and requires the exact final output hash:

```bash
python -m srafeguard.reference_prepare \
  --input-dir /path/to/external/allele_freq \
  --output-dir /path/to/new/external/five_population_build
```

The extraction function is identical to the final manuscript reference producer.
Its original complete four-population reconstruction controls passed before the
validated five-population build. This release reuses those approved inputs and
functions; it does not need private cohort evidence or the original old marker
file. The portability wrapper was tested on an artificial VCF, not rerun over the
full22-chromosome universe in this packaging task. Any final hash mismatch stops
use; inspect formatting/compression/version differences rather than accepting an
approximate panel or substituting another release.

Filters retain the original autosomal/biallelic ACGT-SNV rules, source frequency
format checks and selected-coordinate ambiguity exclusions. Delta AF uses the
maximum pairwise difference across AFR,AMR,EAS,EUR,SAS. Missing AF can be restored
as true zero only with explicit AC=0 and AN>0. Missing/zero denominators remain
unavailable; nonzero AF/AC/AN consistency is checked at the original source float
precision. Merely rescoring an old four-population-filtered table is incomplete.

Table columns: Chromosome, Position (1-based), REF, ALT, AF_AFR, AF_EUR, AF_SAS,
AF_EAS, AF_AMR, Max_pairwise_AF_delta, Most_separated_populations. Last field is
reference diagnostics, not a person/population-origin assignment. The default
profile verifies the final resource SHA. A custom table/profile is not a claim of
frozen manuscript equivalence.

## ClinVar

Frozen resource: CHM13-coordinate lift of ClinVar20220313:
`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz`

Retained reference-hash record:
`a551d3339915c038aaf261200bbe1dc7d25203a44ea6141eac8efc07f9b66d8c`.
This is an old frozen snapshot, not the current ClinVar clinical interpretation.
Download externally and verify the snapshot/checksum; index with tabix or bcftools.

```bash
bcftools index --tbi /path/to/external/chm13v2.0_ClinVar20220313.vcf.gz
```

VCF POS is1-based, BED mask is0-based half-open. Biallelic SNVs are retained;
physical strict-fragment territory determines queried regions, but actual HQ
bases are required for observations. P/LP labels are annotations, not patient
phenotypes or donor genotypes. Do not substitute a live ClinVar release silently.

The bundled four-interval exclusion mask SHA-256 is
`2e1a46b88a2a53908a6966f86a6451f75a9ea4f996f19c49cb08f5e403098b96`.
Its labels are generic recurrent-hotspot method tags, not identifiers or variants.
No individual profile, run linkage or allele counts occur in the BED.

## External software licensing

Fresh scientific environment pins NumPy2.4.6 (BSD-3-Clause), pandas3.0.5 (BSD-3-Clause),
PyYAML6.0.3 (MIT), pytest9.1.1 (MIT), minimap2 2.31 (MIT), samtools/bcftools1.24
(package-reported MIT), SeqKit2.13.0 (MIT) and SRA Toolkit3.4.1 (Public Domain).
Python and transitive binary/library licenses remain applicable; inspect exact
package notices before redistribution. No tool or resource is relicensed by this
software's MIT grant. Exact installed package/build/license metadata is retained
in the separate release review, not bundled executables.
