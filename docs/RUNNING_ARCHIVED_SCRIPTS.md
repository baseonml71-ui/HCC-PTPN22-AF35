# Running archived scripts

Scripts reflect frozen historical analysis implementations developed in sequential stages. They are not presented as one end-to-end executable pipeline. NO_ARCHIVAL is acceptable for archive inclusion. Read release/SCRIPT_ARCHIVE_INDEX.tsv before attempting execution.

Input paths may need adaptation to local downloads, and historical helpers, private administrative records, self-hashes and freeze guards may be unavailable. Obtain original datasets from data/SOURCE_DATA_MANIFEST.tsv. Do not fabricate absent protocols, silently bypass scientific locks, choose whichever wildcard match happens to exist, or treat an adapted execution as the original historical run.

The recoverable broad order is source acquisition and molecular preparation, outcome-blind definition/projection freeze, eligibility/endpoint linkage, patient-level analysis, and manuscript display. Individual modules document more specific ordering in source. No complete global execution order is asserted.

Manuscript-facing source tables are provided for direct result verification. They preserve retained numerical strings, with direct subject/sample/cell/region identifiers removed or replaced by neutral plot-unit labels when pairing is needed. Coordinate and program-vector exports retain native row alignment. Plot labels are not external patient identifiers; no re-identification key is distributed. Rows remain biological/nested units as documented; de-identification does not make cells independent replicates.

Some frozen scripts recursively enumerate a workspace or copy source/protocol records into logs. Run only in a separate directory containing intended public inputs; inspect output destinations before execution so unrelated local files cannot enter output bundles. RC3 did not run these scripts. It neither resolves all historical read expressions nor grants a public redistribution license for every referenced input.

GSE287319 branch was not used for the final manuscript. The three mixed phase15b scripts remain labeled HISTORICAL_SCRIPT. Code-only materialization is an inspection utility, not an approved manuscript analysis entrypoint. The immutable rc1 sanitized sources remain under scripts/frozen_originals; RC2 packaging copies under scripts/current are historical duplicates.
