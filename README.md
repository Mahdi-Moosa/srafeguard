# SRAFEGUARD 0.1.0

SRAFEGUARD quantifies recoverable human-reference-aligned sequence and
genetic-information endpoints in nominally microbial short-read libraries. It
reports conservative human evidence, population-marker recovery, reference-population
affinity and clinically annotated alternate-allele observations.

Version **0.1.0** is an **audit framework**, not a read-depletion tool. Hostile is
an external depletion program; SRAFEGUARD can audit original and externally filtered
libraries separately.

## What it measures

1. Conservative coordinate-deduplicated human evidence and chromosome representation.
2. Clean population-marker loci, five-population likelihoods and block-bootstrap support.
3. High-quality ClinVar-locus observations, observed ALT, candidate P/LP ALT and
   computational placement support, with paired applicability reported explicitly.
4. Per-library machine-readable summaries and optional detailed validation evidence.

These are read-derived observations, not validated diploid genotypes, contributor
identities, disease phenotypes or ancestry assignments. Reference-population affinity
does not establish ethnicity, nationality, geographic origin or identity. Chromosome
representation is not base-by-base genome coverage or a complete human genome.

SRAFEGUARD does not call Hostile, select historical cohorts, build bacterial assemblies,
perform whole-genome donor genotyping, or detect splice junctions. The manuscript's
RNA splice analyses used separate specialized workflows. Genomic alignment audits
can measure genome-alignable signal in RNA reads, but are not splice-aware RNA surveys.

## Installation

Linux is the tested platform. Use the pinned Python 3.14 Conda environment and
install from the repository root:

```bash
conda env create -f environment.yml
conda activate srafeguard-release
python -m pip install .
srafeguard --help
python -B -m pytest -q -p no:cacheprovider tests
```

Dependencies: Python, NumPy, pandas, PyYAML; minimap2, samtools, bcftools, SeqKit,
SRA Toolkit. pytest is for tests. Exact scientific-environment versions are recorded
in environment.yml. External programs must be on PATH: the frozen strict callers
use minimap2/samtools/seqkit directly even where other modules accept tool overrides.
Kraken and bacterial assembly adjudication are not part of this framework's workflow;
unused helper functions in the frozen paired module are retained for source fidelity.
The separate WGS contamination pipeline has a different scope.

Installation from source, wheel and source distribution has been tested in isolated
environments. See [validation](docs/validation.md) for test coverage and limits.
Synthetic tests do not validate every real SRA architecture or reference resource.

## References and configuration

Large references are external, never installed by pip or bundled. Prepare:

- CHM13v2.0 analysis-set FASTA with masked Y PARs and rCRS, plus an `-x sr` index.
- The **complete five-superpopulation-ascertained** CHM13 1KGP marker table.
- The frozen CHM13-coordinate ClinVar20220313 VCF and its tabix index.

The small study-curated recurrent-hotspot BED is included as a reference-only
exclusion mask: four intervals, no alleles, QNAMEs, accessions or person-level records.
Its exact intervals are scientific method criteria, not an individual human profile.
The included `python -m srafeguard.reference_prepare` builds the exact final
marker table from checksum-verified source chromosomes. Its full reference rebuild
was not repeated during software validation. External retrieval, versions and licensing
are in [reference preparation](docs/references.md).

Copy `configs/example.yaml` to a user-owned local configuration. Set external paths;
relative resource paths resolve against the YAML directory. Defaults and resource
manifests are installed package data. There is no implicit research-tree fallback.
`--config` takes precedence over the documented SRAFEGUARD_CONFIG variable; no
current-directory local file is silently loaded.

```bash
srafeguard resources check --config configs/local.yaml
```

The check reports availability; production audit additionally verifies the final
marker-resource checksum in the default FINAL_FIVE_POPULATION profile. Population
thresholds remain configurable. Below Delta AF .10 needs a different complete CUSTOM
resource because the final table starts at .10. The default SRAFEGUARD threshold
is .30, not the .10 primary threshold of manuscript Figure2a. The source-scoring
population order is retained to preserve deterministic tie handling.

See [configuration](docs/configuration.md) for active flat keys, unsupported aliases
and explicit CUSTOM profiles. Tests use a clearly labeled wholly artificial profile.

## Inputs and supported architectures

- Local paired FASTQ/FASTQ.gz: two synchronized short-read files.
- Local single FASTQ/FASTQ.gz: one short-read file; single-read scientific branch.
- Local SRA: layout/platform inferred using READ_TYPE metadata; explicit overrides
  available. No archive read acquisition is performed by SRAFEGUARD.
- Illumina and BGI/MGI paired data use the validated physical-pair caller.
- Ion Torrent single data use the frozen single-read caller. Known long-read
  architectures are rejected. An unknown-platform FASTQ remains the user's
  responsibility to identify as short-read data.

Direct FASTQs receive streaming structural/pair-name validation. Original inputs
are never modified. Local SRA uses vdb-validate, fasterq-dump and the biological
READ_FILTER rule. Only paired biological indices [1,2] are supported in this SRA
path; nonstandard technical/index interleaving is rejected. Use a separately
validated normalization procedure and FASTQ inputs for those cases. This differs
from the newer primary WGS pipeline and is not silently harmonized.

Use the local SRA filename stem as its sample label: the retained extractor expects
matching fasterq output names. Arbitrary aliases are supported for FASTQ audits.
Common FASTQ-name and resource caveats are listed in [limitations](docs/limitations.md).

## Single-library execution

```bash
srafeguard audit --fastq1 /path/to/R1.fastq.gz --fastq2 /path/to/R2.fastq.gz \
  --sample LIBRARY_A --config configs/local.yaml --output workspace/library_a

srafeguard audit --fastq1 /path/to/single.fastq.gz --layout single \
  --platform ION_TORRENT --sample LIBRARY_B --config configs/local.yaml \
  --output workspace/library_b

srafeguard audit --sra /path/to/run.sra --config configs/local.yaml \
  --output workspace/run
```

Sample labels must be simple safe filenames. Strict predicates are fixed: primary
alignments, MAPQ>=30, query coverage>=80%, identity>=98%; paired data require both
mates, proper pairing and the same chromosome, with outer-coordinate deduplication.
Single data use chromosome/start/end deduplication, ignoring strand. Coordinate
deduplication is not proof of independent original molecules. Strict threshold
configuration edits are rejected because the retained callers do not implement
arbitrary strict thresholds. CIGAR/NM details remain those of the frozen callers.

## Batch and re-analysis

Batch TSV columns: `sample`, exactly one of `sra`/`fastq1`, optional `fastq2`,
`reject_source_sra`, `layout`, `platform`. Paths identify user-local files. Never
publish a private manuscript cohort manifest.

```bash
srafeguard batch --manifest /path/to/manifest.tsv --output-root workspace/batch \
  --config configs/local.yaml
srafeguard reanalyze --results workspace/library_a --output workspace/library_a_reanalysis \
  --min-bq 35
```

Re-analysis reads saved effective configuration and original input paths, records
parentage, and writes a new result directory. Original input files must remain
available; retained BAM alone is not its input. `--min-mapq` and `--min-bq` change
population/ClinVar observation QC, not the frozen strict-human caller. Explicit
layout/platform overrides are forwarded to the new analysis.

Batch resume skips completed analysis and cleanup metadata plus privacy_summary.tsv; it does **not**
compare config/reference fingerprints. Use a new output root or reanalyze after
changing configuration, not an assumed cache invalidation. Batch is sequential;
module-level capture/single-read adapters are not thread-safe for concurrent
in-process audits. Independent processes need separate workspaces/sample names.

## PRE/POST use

Audit original and externally filtered files separately with matched labels under
different directories. `--reject-source-sra` can apply the original biological
READ_FILTER policy to externally filtered FASTQs; counts need not match the full
original SRA because external depletion already removed reads. Hostile filtering
and manuscript-specific matched comparison scripts are not included. There is no
`compare` or `hostile` CLI command. Compare matching summary fields locally;
preserve layout and applicability, and do not interpret an unassessed endpoint as zero.

## Output and interpretation

The main `privacy_summary.tsv` uses public report schema `srafeguard-public-2`:

| Layer | Interpretation |
|---|---|
| Strict_units / Counting_unit | absolute deduplicated physical pairs or single reads |
| Human_chromosomes_represented | reference sequences with strict evidence |
| Population_marker_positive | at least one clean usable marker |
| Population_marker_eligible | minimum clean loci met (default20) |
| Population_best_affinity | full-data five-group likelihood winner |
| Population_best_bootstrap | support for that exact full-data winner |
| Population_affinity_recoverable | >=20 loci and winner support>=.90 |
| ClinVar_sites_with_HQ_base_observation | high-quality observations at annotated loci |
| ClinVar_ALT_sites_observed | observed alternate states after hotspot QC |
| ClinVar_candidate_P_LP_ALT_sites | observed ALT with conservative P/LP annotation |
| Pair_supported_P_LP_ALT_sites | computational pair-supported ALT sites, or N/A |
| Single_read_placement_supported_P_LP_ALT_sites | separate single-read locus support, or N/A |

`Strict_human_fragments` is retained as a legacy count name and means single reads
for single-end input; always inspect Counting_unit/Evidence_unit. SRAFEGUARD's
summary is not the primary pipeline's FPM table. It does not add a new denominator
or pretend the two input-normalization paths are interchangeable.

Affinity uses log(AF) for ALT/log(1-AF) for REF, clipping to [.005,.995], with
1,000 resamples of observed chromosome-qualified 1Mb blocks. Missing five-group
frequencies and MIXED loci are excluded; true AF0 is not missing. The bootstrap
supports the full-data winner, not whichever population wins most replicates.
It is a stability assessment of observed blocks, not ancestry or independent
biological sampling. See [methods](docs/methods_audit.md) for exact tie and
observation semantics.

Read-pair-supported P/LP ALT requires expected-locus remapping, retained ALT,
qualifying MAPQ/proper-pair placement and bounded competing-alignment checks.
The mate need not overlap/carry ALT. For single-end data pair support is **N/A**;
its optional single-read placement result remains a distinct endpoint. Neither
is clinical validation or a confirmed donor genotype. Disabled modules are
NOT_EVALUATED/N/A, not measured negative outcomes.

## Privacy and retention

Default: private file permissions, per-library summaries/config/provenance retained;
detailed BAM/FASTQ/QNAME/marker/ClinVar tables and mapping logs removed after audit,
including failure cleanup. Summaries and labels themselves can be sensitive.

```bash
# Explicitly retain sensitive validation evidence:
srafeguard audit --fastq1 /path/to/R1.fq --fastq2 /path/to/R2.fq \
  --sample LIBRARY_A --config configs/local.yaml --output workspace/retained --keep-evidence
```

Aggregate affinity and ClinVar reports are retained by default; per-locus and
read-level evidence requires explicit opt-in. `--no-keep-evidence` explicitly
requests cleanup, and an omitted flag uses the configured policy. Temporary scratch
is removed even when final evidence retention is requested. Cleanup failures are
reported separately from analytical completion; original inputs are protected.
`runtime.keep_intermediates=True` is unsupported; use the retention flags instead.

For audit and batch, explicit CLI retention flags override the supplied configuration.
For re-analysis, precedence is: explicit retention flag, current `--set` values,
current `--config` overlay, inherited effective configuration, then defaults.
The effective configuration and metadata record the policy actually applied.
Resume cannot restore evidence deleted previously; re-analysis requires the original
inputs. See [reporting and retention](docs/reporting_compatibility.md) and
[privacy](docs/privacy.md).

No original cohort roster, selection seed, accession mapping or human profile is
included. The fixed bootstrap seed is a resampling setting only, not a historical
cohort-selection seed. The program does not select cohorts. Pseudonyms, summary
reports and passing disclosure scans do not guarantee anonymity.

## Executable safe example

```bash
python examples/synthetic_workflow.py --workdir /tmp/srafeguard-artificial-example
```

Choose a new directory outside this source tree. It generates wholly artificial
references, reads, AF values and VCF annotation labels; uses actual minimap2,
samtools and bcftools; exercises paired/single audit, retention, batch/resume and
reanalyze. It is not an actual human allele profile or clinical validation.

## Reporting compatibility

The default `srafeguard-public-2` schema has 44 fields. Exact 35-field compatibility
with historical manuscript analysis is available through `reporting.schema: legacy-1`;
the public companion preserves explicit applicability and evaluation status. Legacy
fields containing `validated` describe historical computational endpoints, not
clinical validation. See [reporting compatibility](docs/reporting_compatibility.md)
for configuration, nullable parsing, supporting-read/mate counters and missing-value
semantics.

## Resource requirements and limitations

Reference identities, checksums and external data terms are documented in
[reference preparation](docs/references.md). Large human indexes can require roughly
10 GB or more RAM; the full population panel and temporary BAMs add memory and disk
requirements that depend on library size. Start with small thread counts and separate
user-owned scratch. No fixed throughput or universal RAM guarantee is made.

Validation includes real-tool synthetic workflows and saved scientific-output
regression checks. It does not constitute complete historical cohort reprocessing,
full real-SRA integration or clinical validation. See [validation](docs/validation.md)
and [limitations](docs/limitations.md). Abrupt process termination or power loss can
prevent cleanup; inspect private workspaces before reuse.

## License and citation

Original SRAFEGUARD code is licensed under the [MIT License](LICENSE),
Copyright (c) 2026 Mahdi Moosa. External tools, human reference genomes, 1KGP and
ClinVar data retain their own licenses and provider terms; see
[references](docs/references.md).

Please cite SRAFEGUARD, its author and version using [CITATION.cff](CITATION.cff).
