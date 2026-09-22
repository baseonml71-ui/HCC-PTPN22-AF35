# HCC PTPN22–AF35 analysis — local RC2

**Release gate: NOT_READY.** Stop state: HCC_PUBLIC_REPOSITORY_V1_RC2_AUTHOR_REVIEW. This repository has not been published, pushed or uploaded.

This candidate preserves 59 sanitized frozen rc1 source files, 59 current source candidates, fixed AF35 membership, historical software records, 43 aggregate publication tables and a source-result evidence index. It has **zero approved current analysis entrypoints**. Removal of eight administrative private-file copies did not close scientific input contracts. Public raw datasets must be obtained from their original repositories; accession records alone do not establish complete reproduction instructions.

The [release gate](PUBLIC_REPOSITORY_RELEASE_GATE.md), [traceability audit](release/PUBLIC_RESULT_TRACEABILITY.tsv), [source map](release/PUBLIC_ANALYSIS_SOURCE_MAP.tsv) and [unresolved actions](release/UNRESOLVED_AUTHOR_ACTIONS.tsv) are authoritative for this candidate. No full one-command or source-to-result reproducibility claim is made.

## Scientific scope

The manuscript examines a PTPN22-associated activation-feedback T-cell state across patient-level clonotype, spatial and clinical contexts. Patient/sample is the biological replicate; cells and regions are nested observations. Association does not establish causal PTPN22 control of AF35, a universal ICI predictor or a clinical cutoff. See [scientific limitations](docs/LIMITATIONS.md).

## Reproducibility levels

**Level 1 — runnable from included derived data:** package integrity checks and inspection of 43 exact aggregate source tables are available. No scientific analysis pipeline has been demonstrated runnable from those tables. AF35 membership is directly verifiable from the included 35-gene vector.

**Level 2 — requires public source datasets:** intended analytical reconstruction uses the [source manifest](data/SOURCE_DATA_MANIFEST.tsv) and [acquisition guidance](data/DOWNLOAD_INSTRUCTIONS.md). No analysis currently qualifies as a verified Level 2 pipeline because deterministic module contracts remain incomplete.

**Level 3 — source provenance with raw data absent:** source identities and hashes are recorded where verified, but most chains remain TRACE_INCOMPLETE. They must not be read as completed Level 3 reproducibility. Original raw and individual-level data are not redistributed. Exact availability and upstream producer uncertainty are retained in the [derived-input ledger](release/DERIVED_INPUT_LEDGER.tsv).

## Verification and review order

```sh
python scripts/utilities/verify_release.py
python scripts/utilities/inspect_prerequisites.py
```

The first command checks payload integrity and AF35 membership only. The second intentionally exits with code 2 because reproduction blockers remain. Neither runs analyses or downloads data. The optional scripts/utilities/materialize_workspace.py copies candidates into the original module layout under ignored data/local_workspace for code inspection only.

Review family documentation in order:

1. [Anchor](scripts/current/01_anchor/README.md)
2. [State](scripts/current/02_state/README.md)
3. [AF35 measurement](scripts/current/03_AF35_measurement/README.md)
4. [Clonotype](scripts/current/04_clonotype/README.md)
5. [Spatial](scripts/current/05_spatial/README.md)
6. [Clinical](scripts/current/06_clinical/README.md)
7. [Perturbation](scripts/current/07_perturbation/README.md)

These are review entrypoints, not analysis launchers. The [packaging boundary](docs/RC2_PACKAGING_BOUNDARY.md) explains frozen-original identity, exclusions and static audit limitations. No missing scientific protocol was reconstructed from memory.

## Definition, provenance and data

[AF35 membership](config/AF35_genes.txt) remains 35 ordered unique genes, unsigned and equally weighted, excluding PTPN22. See [definition](docs/AF35_DEFINITION.md) and [origin](docs/AF35_PROVENANCE.md). The original exact vector and creation patch were recovered; a gene-by-gene reduction rule was not recovered. No definition, statistic, estimand, manuscript or figure was changed during RC2.

The [table bindings](derived_data/manifests/SOURCE_TABLE_BINDINGS.tsv) and [frozen-to-current map](release/FROZEN_TO_PUBLIC_SCRIPT_MAP.tsv) preserve hashes. Figures, decorative assets and individual-level inputs remain excluded. Original source terms still apply. The [environment records](environment/README.md) are incomplete historical records; the topology SciPy version is unknown.

## Citation and license

[CITATION.cff](CITATION.cff) records the supplied six authors and manuscript title. [Zenodo metadata](ZENODO_METADATA_DRAFT.md) remains an unpublished draft. No repository URL or release DOI has been invented. [MIT](LICENSE) covers original code and accompanying documentation, not third-party datasets or assets. See [exclusions](release/EXCLUDED_FILES.tsv).
