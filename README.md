# HCC PTPN22–AF35 analysis

## Overview

Code and scientific provenance for a PTPN22-associated activation-feedback T-cell state in hepatocellular carcinoma (HCC). Planned repository: HCC-PTPN22-AF35; software version 1.0.0. This is an unpublished local candidate. See [release readiness](PUBLIC_REPOSITORY_RELEASE_GATE.md).

## Scientific question

Can a fixed distributed T-cell-state readout connect patient-resolved clonotype context, conditional tissue organization and pretreatment clinical association?

## Main findings

PTPN22 and AF35 mark expanded, persistent and baseline later-observed clonotype contexts. State adjustment substantially attenuates these contrasts. Spatial and pretreatment clinical analyses extend the state framework with cohort-dependent evidence. Perturbations affect selected activation dimensions without uniformly controlling AF35. No causal PTPN22-to-AF35 relationship, ICI-specific predictor or clinical cutoff is established.

## Repository scope

This repository provides analysis code, the fixed AF35 definition, environment records, provenance documentation and selected aggregate derived source tables used to support the reported analyses. Original public datasets should be obtained from their source repositories using the accession information provided. It is not a demonstrated fully reproducible raw-data pipeline.

## Analysis architecture

Anchor selection → patient-level state architecture → fixed AF35 measurement → clonotype/temporal contexts → conditional spatial organization → clinical association; orthogonal perturbations define mechanistic limits. See [analysis order](docs/ANALYSIS_ORDER.md) and [entrypoints](provenance/CURRENT_ANALYSIS_ENTRYPOINTS.tsv). Patients or independent biological replicates are the inferential units.

## AF35 definition

[AF35_genes.txt](config/AF35_genes.txt) contains 35 ordered unique genes, fixed and unsigned with equal weights, excluding PTPN22. Platform-specific score implementations are documented in [AF35 definition](docs/AF35_DEFINITION.md). The exact original specification and creation patch were recovered; the gene-by-gene reduction rule was not retained. GSE238264 is specification context and Wu the independent pretreatment spatial extension. See [public provenance](docs/AF35_PROVENANCE.md).

## Dataset sources

The [source manifest](data/SOURCE_DATA_MANIFEST.tsv) identifies accessions, original studies, roles, overlaps and terms. [Download instructions](data/DOWNLOAD_INSTRUCTIONS.md) point to GEO, ENA and the original Wu Zenodo record. No external raw or individual-level data are redistributed.

## Reproducing analyses

Start with a static package check using Python 3 (standard library only):

```sh
python scripts/utilities/verify_release.py
python scripts/utilities/inspect_prerequisites.py
```

The first command checks package integrity. The second reports known execution gaps and intentionally exits with code 2 while public reproduction prerequisites are unresolved. Neither command runs biological analyses or downloads data.

For inspection of the original module layout, optionally run:

```sh
python scripts/utilities/materialize_workspace.py
```

This only copies code to ignored data/local_workspace. It does not make unavailable inputs available. Do not run archived analyses until the module's source data, locked inputs, environment and public protocol dependencies have been resolved. Follow [reproducibility notes](docs/REPRODUCIBILITY_NOTES.md); do not bypass historical hashes or response-blind locks.

## Directory structure

docs contains scientific scope and methods boundaries; config contains fixed definitions; scripts contains categorized analysis sources and utilities; environment contains historical version records; data contains accessions and acquisition guidance; derived_data contains selected aggregate results and bindings; provenance contains lineage and entrypoint maps; figures documents asset exclusions; supplement indexes source tables; release contains notes, exclusions, manifest and checksums.

## Software environments

See [R records](environment/R_environment.txt), [Python records](environment/Python_environment.txt) and module-specific session_info. The historical SciPy version for the topology module was not recoverable. No complete environment lockfile is fabricated.

## Derived outputs

The selected aggregate tables in derived_data/publication_source_tables are byte-identical to frozen sources. [Bindings](derived_data/manifests/SOURCE_TABLE_BINDINGS.tsv) preserve the checksums. Individual-level parts are omitted; [exclusions](release/EXCLUDED_FILES.tsv) explain the selection. No figure binaries are included while decorative-asset rights remain unresolved.

## Data availability

Source datasets remain with their original repositories and terms. This candidate is local only; no public repository URL, release DOI or archive DOI exists yet. Accessions refer to source datasets, not a DOI for this software.

## Citation

Use [CITATION.cff](CITATION.cff) for software authorship and the associated manuscript title. Repository URL and DOI are pending publication and are intentionally absent from identifier fields. See [Zenodo metadata draft](ZENODO_METADATA_DRAFT.md).

## License

[MIT](LICENSE) applies to original repository code and its accompanying documentation. It does not change the licenses or terms of external datasets, third-party publications, third-party software or externally sourced materials. Inclusion of aggregate study summaries does not authorize redistribution of their original individual-level inputs.

## Known limitations

See [scientific limitations](docs/LIMITATIONS.md) and [reproducibility notes](docs/REPRODUCIBILITY_NOTES.md). The source archive retains unavailable historical protocol dependencies and some mixed historical/current branches. Static file verification is not a successful biological rerun. The public release gate is NOT_READY until the documented input contracts and source lineage are resolved.
