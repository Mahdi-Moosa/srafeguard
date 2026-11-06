# Scope and unresolved validation limits

- Full CHM13/ClinVar/1KGP clinical integration was not freshly rerun here. Saved
  frozen-score receipts were checked; actual installed tools also ran end-to-end
  on wholly artificial references, variants and reads. These are different levels
  of evidence, not a manufactured clinical/ancestry validation.
- The portable marker builder retains the frozen extraction function and22 source
  hashes, but a full8.83GiB rebuild was not repeated. It stops on output hash mismatch.
- Local SRA tools are installed, and safe mocked input/REJECT tests exist; actual
  local-SRA end-to-end integration was not performed. Paired indices beyond[1,2]
  are explicitly unsupported, unlike the newer primary WGS pipeline. Use normalized
  FASTQs rather than silently selecting biological/technical files incorrectly.
- Use the SRA filename stem as sample name. The original SRA extractor derives
  expected fasterq filenames from that label. Arbitrary FASTQ aliases are supported.
- Direct paired FASTQ names must have identical first whitespace-delimited tokens.
  /1,/2 suffixes are rejected before audit because the retained SeqKit export
  assumes mapper/FASTQ identifiers agree. Normalize suffixes consistently outside
  this tool; no new read-naming transformation was silently introduced.
- Batch resume tests analytical and cleanup completion plus summary existence, not configuration or
  reference equivalence. Use new output roots or reanalysis after changes.
- Capture and single-read adapters temporarily patch module functions; concurrent
  in-process library audits are unsupported. Batch is serial. Separate processes
  must use distinct scratch/workspaces and sample labels.
- Runtime outputs remain sensitive even when details are cleaned. Summaries,
  effective configs and provenance can identify inputs and recoverable information.
  External shell logs are the user's responsibility; set `umask 077` before capture.
- Missing resources fail the audit. Disabled modules are NOT_EVALUATED/N/A. Pair
  support on single-end data is N/A, while independent-read placement remains a
  separate scientific result. Empty support is not clinical absence.
- No Hostile, bacterial assembly, RNA splice, diploid genotyping or identification
  function is implemented. PRE/POST use means matched separate audits, not an
  automatic depletion/comparison mode. No historical cohort can be reconstructed
  from withheld sampling seeds or rosters in this distribution.
