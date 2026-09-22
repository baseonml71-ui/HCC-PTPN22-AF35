# Manuscript-facing source tables

Final manuscript values can be verified against the included manuscript-facing source tables. Table exports preserve retained numerical strings and remove/neutralize direct identifiers. No re-identification key is included. Original and derivative hashes and removed columns are in [bindings](../../release/FINAL_SOURCE_TABLE_BINDINGS.tsv).

The [main table index](../../release/MAIN_TABLE_INDEX.tsv) maps every included main source table to Figure, Panel, Description, Units, Biological replicate and Source dataset. Supplement sources are in the [supplement index](../../release/SUPPLEMENT_PROVENANCE_INDEX.tsv). Units vary by statistic and retain their source column definitions; cohort-native scores are not interchangeable.

GSE206325 coordinate, marker and UCell tables retain identical row order after cell IDs are removed. GSE235863 atlas/program tables retain their frozen native row order. Shared neutral plot units preserve the seven-patient longitudinal pairing and the GSE302495 endpoint/score join. Do not link plot units across unrelated cohorts or treat them as source patient identities.

F4B_COORDINATE_SOURCE_RECORD is a public-source pointer plus verified display count/hash, not redistributed CODEX coordinates. Original source data remain at the linked Zenodo record. H&E and figure binaries are absent.
