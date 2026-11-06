# Privacy and public-report policy

The release includes reviewed source, synthetic tests and a reference-only four-
interval hotspot BED. No real allele vectors, variant observations, human reads,
private labels/crosswalks/selection seeds or original Git metadata are distributed.
Default file creation uses umask077. Direct FASTQs receive streaming structure and
pair-name checks; original input files are never deleted or rewritten.

The research program default removed major read evidence but retained per-locus
marker and ClinVar/QNAME tables. The public candidate adds a separate explicit
cleanup layer: with keep_evidence false it deletes only explicitly enumerated
owned sensitive detail files and raw mapping logs, including on failure. Aggregate
population scores, clinical endpoint counts and chromosome summaries are retained.
Unknown files in evidence/log directories are preserved, not recursively deleted.
Keep-evidence opt-in retains
them for controlled validation. The program still leaves per-library summary,
effective configuration and provenance; these can expose sensitive information
and private input paths and must not be committed or publicly uploaded blindly.

Canonical privacy_summary.tsv adds explicit applicability and replaces ambiguous
validated aliases with pair-supported versus single-read-placement-supported
fields. Internal historical variable names/classifiers are unchanged. Single-end
pair support and disabled modules are N/A, never a substituted biological zero.
A schema label distinguishes this first public report from historical private
exports. Optional retained legacy analytical files retain their older column names.

Reference-population affinity and clinical annotation do not establish identity,
donor genotype, ancestry or disease status. Numerical fingerprints, sample labels,
input filenames and cross-table joins can create linkage risks. Pseudonymization
or an aggregate summary is not a guarantee of anonymity. Secure storage, access
control and disclosure review remain the local user's responsibility.

The fixed population bootstrap seed is not an original cohort sampling seed and
cannot reconstruct withheld cohorts. No cohort-selection module exists here.
No secrets, credentials, external .git objects, actual human genetic profiles or
private selected lists may enter future source archives, wheels or Git staging.


The authoritative retention/output policy is now [reporting compatibility](reporting_compatibility.md): aggregate affinity/ClinVar files are retained, and only explicit owned detail files are removed. Unknown files in evidence/log directories are preserved.
