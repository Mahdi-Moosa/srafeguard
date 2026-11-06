# Reporting and retention compatibility

The default report schema is `srafeguard-public-2` (software version remains 0.1.0). Scientific classifiers and calculations are unchanged. All reports are local potentially sensitive outputs, not anonymous publication artifacts.

## Retention

Without `--keep-evidence`, the program removes only explicitly enumerated owned detail files: strict BAM/index, strict FASTQs/QNAMEs/fragment-coordinate tables, marker observations, raw/filtered per-locus ClinVar tables, read-placement evidence and raw mapping logs. Unknown files are never recursively deleted. Only empty evidence/log directories are removed. The exact allowlist is `release_policy.SENSITIVE_FILES`. Original input paths are protected.

Default retained aggregates include `privacy_summary.tsv`, `population_affinity.tsv` (all likelihood and bootstrap fractions), `clinvar_exposure.tsv` (endpoint/site/base/support counts), `single_read_clinvar_validation.json` where produced, `evidence/chromosome_counts.tsv`, contamination summaries, effective configuration and provenance. These files contain per-library aggregate measurements, not locus coordinates, QNAMEs, allele vectors or variant identities. They can still permit linkage; handle responsibly and do not automatically publish them.

Use `--keep-evidence` for deliberately retained human evidence. This retains owned detail and raw logs, not all temporary scratch alignments. `runtime.keep_intermediates` remains unsupported. Completed metadata separates `analytical_status`, `cleanup_status`, `cleanup_errors`, `evidence_retention_requested`, `evidence_retention_status` and `evidence_retained`. Runs are completed only after successful cleanup; errors fail closed. Failure `diagnostic_summary.json` records phase/error type without raw command, path, QNAME or sequence text. Private run metadata/configuration retain input/resource paths. Failure messages are sanitized by default; raw exception detail is retained only with explicit keep-evidence.

Batch resume requires completed analysis and cleanup metadata. It does not compare configuration/resource fingerprints or restore evidence deleted previously. Requesting keep-evidence on a skipped completed run cannot recover files: rerun from original inputs into a new output directory, or use `reanalyze`. Re-analysis reruns original input paths from parent provenance; aggregate reports or retained BAMs alone do not replace original inputs. Preserve those inputs separately.

## Public schema

The 41-field public-1 report gains three directly measured fields (44 total): `Pair_supported_P_LP_supporting_mates`, `Single_read_placement_supporting_reads`, `Single_read_placement_evaluated`. No counts are inferred from site counts. Existing `Strict_units`, `Counting_unit`, marker-positive/eligible/status, pair applicability/evaluation and paired/single supported-site fields remain. Pair support for single-end input is literal `N/A`; independent-read placement has its own fields. Zero means an evaluated negative. Disabled module endpoints and disabled confirmation are `N/A`, not zero. Known MAPQ/BQ/chromosome-scope/mask configuration values are preserved even when screening is disabled. `ClinVar_evaluation_status=EVALUATED` refers to screening; check placement evaluation flags separately.

The old four fields containing `validated` are not scientific clinical-validation claims and are absent from the default public summary. Their exact original information is available through explicit legacy output.

## Exact historical compatibility

Add to an otherwise valid configuration:

```yaml
reporting:
  schema: legacy-1
  legacy_export: true
```

This writes the original 35 fields and original values into `privacy_summary.tsv` (also `privacy_summary_legacy.tsv`), and the improved report into `privacy_summary_public.tsv`. Alternatively use `schema: public-2` with `legacy_export: true` to keep the public primary report and request only the legacy sidecar. No remapping or scientific recalculation differs between schemas.

Legacy fields preserved exactly: `ClinVar_validated_P_LP_ALT_sites`, `ClinVar_validated_P_LP_supporting_reads`, `ClinVar_validated_P_LP_supporting_mates`, `P_LP_validated_exposure_detected`. Supporting counts are the original summed per-site counts and may count the same read at multiple sites; they are not unique donors, molecules or necessarily unique reads. The historical single adapter aliases independently supported reads into the mates counter; its value must not be interpreted as physical pair confirmation. Paired supporting_reads can be blank as originally produced. Legacy disabled modules preserve historical zero/false/default status behavior, which is ambiguous: always consult the public companion and effective configuration; legacy zeros are not evidence an assay was performed. Legacy fields do not establish clinical validation or donor genotype.

Figure 3/4 historical preparers require legacy-named fields. Run in `legacy-1` mode for those workflows; no keep-evidence is required for their aggregate counts. Five-population migration extraction requires marker observations and/or read evidence, so explicitly request `--keep-evidence` as well. Its historical script is tied to private resources and original analysis settings; preserving files is not a claim that it runs unmodified on arbitrary new inputs.

## Nullable parsing

```python
import pandas as pd
x = pd.read_csv('privacy_summary.tsv', sep='\t',
                keep_default_na=False, na_values=['N/A', ''])
x['Pair_supported_P_LP_ALT_sites'] = pd.to_numeric(
    x['Pair_supported_P_LP_ALT_sites']).astype('Int64')
def nullable_bool(value):
    if pd.isna(value):
        return pd.NA
    return {'true': True, 'false': False}[str(value).lower()]
x['Population_marker_eligible'] = x['Population_marker_eligible'].map(
    nullable_bool).astype('boolean')
```

Text reference-population fields may be empty when no score exists. Never apply `.astype(bool)` to arbitrary strings or missing values. TSV has no enforced embedded types. Consumers should check `Report_schema`; legacy schema has no such column to preserve exact original layout.

## Retention precedence and owned temporary evidence

For audit/batch, explicit --keep-evidence or --no-keep-evidence overrides the supplied config, then installed defaults. For reanalyze: explicit retention flag > current --set values > current --config overlay > inherited effective config > program defaults. Current config resource/temp paths resolve relative to that file. No override inherits the parent's effective policy. The child saved config, effective metadata and actual retention all use the resolved policy; parent evidence is not copied. Resume cannot restore deleted evidence; rerun original inputs into a new output directory.

Owned temporary evidence is registered in owned_sensitive_intermediates.json before writing sensitive contents. This includes the selected-unsorted BAM, isolated samtools sort spill directory, captured original mapping BAM and uniquely created input/mapping/clinical scratch directories. Only those exact registered temporary owners are cleaned; final evidence remains subject to keep-evidence. Temporary scratch is removed even when final evidence retention is requested. Cleanup errors are explicit; aggregate reports remain available. The registry/configuration are private local provenance with paths, not publication exports. Unknown user directories/files and original inputs are protected.
