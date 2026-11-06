# Scientific fidelity and implementation boundaries

SRAFEGUARD software version0.1.0 is verified from the original pyproject; the
scientific profile is named separately. The active CLI is cli -> runner -> pipeline
with local input normalization, contamination, original-alignment capture,
population and ClinVar modules. Older reference/ backups, legacy Human Leakage
Audit README, historical batch manifests and manuscript-specific scripts are not
part of the public entry points.

| Endpoint | Supported implementation |
|---|---|
| Strict paired evidence | original frozen paired caller, both primary proper mates, same chromosome, MAPQ30, coverage80%, identity98%, outer-span dedup |
| Strict single evidence | frozen single caller, MAPQ30/80%/98%, chromosome/start/end dedup ignoring strand |
| Unit and breadth | absolute physical-pair/single-read counts and represented reference sequences, not DNA molecule number or genomic coverage |
| SNP extraction | retained original strict-human BAM, manual CIGAR/base projection, MAPQ30/BQ30, REF/ALT counts collapsed by locus |
| Clean markers | REF-only or ALT-only among usable REF/ALT bases; MIXED excluded; non-REF/ALT bases ignored by the retained extraction function |
| Population model | five-group per-locus log(AF) or log(1-AF), clipped [.005,.995]; missing frequency rows excluded |
| Blocks | chromosome-qualified (Position-1)//1000000 observed blocks; within-library replacement draws preserving block count |
| Support |1000 NumPy default_rng resamples, support for full-data LL winner, >=20 clean loci and >=.90 default recoverability |
| ClinVar screening | physical strict territory, biallelic SNVs, raw and post-hotspot observations, conservative P/LP annotation |
| Pair support | ALT-carrying mate recovers target ALT, primary/proper paired placement, MAPQ30, top score, no qualifying equal-or-better distinct-locus competitor |
| Competitor search | original -N50/-p.5/secondary=yes and500bp distinct-locus window; bounded search, not proof of unique genomic origin |
| Single support | original independent-read placement adapter; not a physical pair and reported separately |

Full-data ties use the original np.argsort(...)[::-1] order; bootstrap ties use
np.argmax. Population order is preserved, so tie behavior can make likelihood
winner and bootstrap-leading group disagree. The canonical support is indexed by
the full-data winner. It is never replaced by maximum bootstrap support.

The fixed seed is a resampling configuration only. Blocks are observed evidence
units, not demonstrated independent original donors/fragments. No genotype,
ancestry, identity, nationality, donor diagnosis or causal source is established.
ALT-supporting pairs need not have the ALT on both mates.

The main WGS pipeline and this framework share strict scientific predicates but
are not interchangeable execution/counting systems. The newer WGS acquisition
supports additional biological read-index layouts. SRAFEGUARD supports local SRA
paired indices[1,2] and single layouts, with its own READ_FILTER normalization;
direct FASTQ uses supplied reads plus optional reject-source normalization.
SRAFEGUARD's absolute-count report does not reproduce a primary-pipeline FPM
without a separately justified input denominator.

## Explicit public-release changes

- Installed defaults/resource manifest, relative path resolution and user cache
  replace development-tree coupling. Source scoring/calling classifiers preserved.
- Frozen strict QC edits and inactive nested config aliases are rejected, rather
  than pretending unsupported settings take effect. Genomic assay only; known
  unvalidated long-read architectures rejected.
- Nonstandard paired SRA read indices now fail early instead of possibly choosing
  technical FASTQ files. No new biological extraction algorithm was invented.
- Provenance stores categorical layout/platform without path conversion.
  Re-analysis forwards saved explicit layout/platform; the original adapter omitted
  these. Tool-version probing ignores failed --version invocations before trying
  alternatives. These are provenance/routing fixes, not scientific caller changes.
- Canonical privacy summary replaces ambiguous validated aliases with explicit
  pair-versus-single applicability and preserves internal legacy scientific fields.
- Default cleanup now removes detailed human marker/ClinVar tables and logs, not
  just reads. Opt-in retention preserves them. Retention-on/off numerical summaries
  match in the artificial end-to-end test.

Original original-alignment population code matches the final migration producer.
Saved private regression receipts and393 score-file checksums support the known
WGS91/39 PRE, POST7/0/0, median68/0, and RNA9/0 results. These are saved scientific
outputs, not outcomes of fabricated fixtures or a fresh historical rerun.

No manuscript or supplementary text was edited. The current package presents an
audit, not depletion. Single-end physical pair support is N/A. Clinical placement
support is not clinical validation. Original historical bootstrap-producer and
other manuscript reconstruction limits are not claimed solved by packaging.
