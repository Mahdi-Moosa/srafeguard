# Active configuration and compatibility

Installed `srafeguard/data/default.yaml` is the effective default. Explicit --config
or SRAFEGUARD_CONFIG supplies overrides; resource paths are relative to the override
file. resources check is an availability report, not reference identity proof.
The default FINAL_FIVE_POPULATION profile pins the complete final marker table and
requires the known human FASTA checksum and the frozen ClinVar checksum when
clinical auditing is enabled. The minimap2 index/FASTA construction lineage and
the VCF/index pairing must also be verified externally. CUSTOM_SYNTHETIC is used only
for fabricated integration assets; an arbitrary custom resource is not a claim of
manuscript concordance.

Active population keys: enabled, extraction_mode=original_alignment, populations,
min_pairwise_af_delta, af_floor, min_mapq, min_bq, block_size_bp, bootstraps, seed,
informative_min_sites, informative_min_bootstrap. Keep all five groups. Defaults
are .30/.005/30/30/1000000/1000/fixed bootstrap seed/20/.90. The original order
AFR,EUR,SAS,EAS,AMR is retained for exact deterministic tie behavior. Membership of
the group set, not alphabetical order, defines five-population scoring.

Inactive duplicate nested population.bootstrap/remap/interpretation and
superpopulations aliases from the research YAML are rejected rather than silently
ignored. Use the active flat keys. Sub-.10 thresholds require a different complete
CUSTOM resource because the final panel has an ascertainment floor of .10.

ClinVar keys: enabled, min_mapq, min_bq, biallelic_snv_only, chromosome_scope,
apply_hotspot_exclusion_mask, validation.enabled/threads/secondary_n/secondary_p/
min_remap_mapq/distinct_locus_window_bp/min_passing_mates. Defaults preserve the
validated classifier, including 50/.5 competing-alignment search, MAPQ30 and a
500bp distinct-locus window. These bounded checks are not a uniqueness proof.

runtime.temp_dir is also passed as runtime.tempdir for the retained modules; the
scientific behavior is unchanged. outputs.keep_evidence or --keep-evidence controls
post-analysis retention. output.retain_* toggles in the old YAML were not wired to
cleanup and are not advertised. runtime.keep_intermediates is disallowed by the
release policy. Fixed strict-human thresholds cannot be altered by configuration;
changes are rejected instead of falsely appearing to take effect.

The public recorder keeps layout/platform as categorical strings, not absolute
paths (the original recorder treated every input value as a path). Re-analysis
copies recorded layout/platform, unlike the original adapter which
omitted those overrides. It still requires the original read files. Batch resume
only checks completed status/summary presence, not parameter/resource equivalence.
