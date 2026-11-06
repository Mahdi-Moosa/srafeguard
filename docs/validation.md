# Validation

This record summarizes bounded validation of SRAFEGUARD 0.1.0. The final suite
passed **77 tests** with no failures, errors or skipped tests: 54 existing tests
and 23 additional cleanup, failure-propagation and retention-precedence tests.
These checks are distinct from historical scientific cohort reprocessing.

## Installation and external tools

Source, wheel and source-distribution installation were tested in isolated Python
and Conda environments. The documented environment uses Python 3.14 and the tool
versions in [environment.yml](../environment.yml). CLI entry points, configuration
loading, reference-path errors, dependency resolution and research-directory
independence were checked. Final wheels and source distributions were also tested
in a fresh joint Python environment with the WGS contamination pipeline.

Use the documented environment: package metadata declares Python >=3.10, but the
current dependency family requires a newer interpreter (>=3.11). Older interpreter
support was not established. Installation checks do not change scientific callers.

## Synthetic workflow and reporting

The [synthetic example](../examples/synthetic_workflow.py) runs actual minimap2,
samtools and bcftools with wholly artificial sequences, reference frequencies and
VCF labels. It exercises paired and single FASTQ audits, machine-readable outputs,
optional evidence retention, batch/resume and re-analysis. Artificial inputs are
not actual human allele profiles or evidence of clinical/ancestry accuracy.

Tests cover strict QC boundaries, chromosome summaries, clean marker eligibility,
five-population likelihoods and bootstrap reproducibility. Support is evaluated
for the full-data likelihood winner. ClinVar screening, candidate/placement-supported
endpoints, single-end N/A, disabled-module status and known configuration values
remain distinct.

The default 44-field public schema and exact historical 35-field export are tested.
Original Figure 3 count/boolean consumers and Figure 4 endpoint expressions were
executed against legacy outputs. Entire figure-generation and RNA preparation jobs
were not repeated. Detailed affinity migration requires explicitly retained
read/per-locus evidence and its original resources; aggregate reports alone do not
satisfy those inputs. See [reporting compatibility](reporting_compatibility.md).

## Retention and technical failures

Earlier regression checks established aggregate-report retention and explicit
legacy exports. Subsequent checks on 5–6 November 2026 covered failed sorting,
partial BAM output, subprocess exceptions, cleanup deletion failure, protected
inputs, effective-policy metadata and re-analysis overrides. Real samtools sort
failure with partial output was tested, rather than relying only on mocked cleanup.
The original independent failure reproductions were checked after correction.

Owned temporary evidence is registered before writing. Successful cleanup verifies
its removal; cleanup failure remains explicit and cannot become analytical success.
Completed aggregate affinity/ClinVar reports remain available. Explicit keep/cleanup
choices and configuration inheritance were tested across audit, batch and re-analysis.
Resume does not restore previously deleted evidence. SIGKILL or power loss can bypass
normal cleanup and requires workspace inspection.

## Saved scientific regression

All **393 saved score files** were rehashed and reconciled with retained validation
records. Frozen aggregate checks include WGS PRE 91/100 marker-eligible and 39/100
affinity-recoverable libraries; POST 7/100 marker-positive, 0/100 eligible and
0/100 recoverable; median clean marker counts 68 to 0; RNA-seq affinity recovery
9 PRE to 0 POST among 18 libraries. This is saved-output verification, not fresh
mapping, extraction or historical bootstrap computation.

Scientific population likelihood/bootstrap code, strict callers, legacy calculations
and clinical classifiers were preserved. Technical corrections concerned subprocess
failure handling, temporary ownership, cleanup and configuration precedence.

## Coverage limits

- Full historical cohorts were not reprocessed.
- Real local-SRA end-to-end coverage is limited; nonstandard paired biological indices
  beyond [1,2] are unsupported. FASTQ validation has separate tests.
- Complete CHM13/ClinVar real-reference integration and the full population-reference
  rebuild were not repeated during packaging. The reference builder retains its source
  checksums and fails on output mismatch.
- The inactive `release.report_schema` label remains `public-1`; operative
  `reporting.schema` and generated default reports use `public-2`.
- Synthetic tests do not demonstrate general sensitivity/specificity, donor identity,
  clinical interpretation or universal robustness to changing tools/resources.

Source, wheel and source-distribution disclosure checks and manual inventory reviews
were performed separately from algorithm tests. No real human read/allele evidence
is distributed as a fixture. Passing scans and aggregate reports do not guarantee
anonymity. See [privacy](privacy.md), [limitations](limitations.md) and
[reference preparation](references.md).
