# Reproducibility notes

## What has been verified

The fixed ordered vector matches the recovered specification. Included aggregate tables are byte-identical to frozen publication sources. Copied Python code parses, and its non-string syntax tree is identical to the source: only path strings were sanitized. No biological code was executed, no cutoff was optimized and no frozen results were regenerated. The release checks verify package files and metadata, not scientific replication.

## Workspace layout

Scripts retain their module-relative layout under scripts/<analysis area>/<original relative location>. The optional utilities/materialize_workspace.py helper copies those scripts into an ignored data/local_workspace directory using provenance/SCRIPT_SOURCE_BINDINGS.tsv. It does not download data or run analysis. Original relative project roots are then resolved from that isolated workspace. Do not execute the archived modules directly from the categorized script folders because sibling-module relationships would differ.

## Inputs and outputs

provenance/SCRIPT_FILE_REFERENCES.tsv lists statically extracted file literals; it is not a complete dependency graph. provenance/CURRENT_ANALYSIS_ENTRYPOINTS.tsv identifies the script, language, environment, and script-defined output context. Dynamic expressions, inherited helper roots and external freeze manifests require module-level review. Expected output directories are created by the archived modules in an isolated workspace; they are not prepopulated with sensitive data.

The public accessions are listed in data/SOURCE_DATA_MANIFEST.tsv. Acquire original inputs under their source terms. Individual-level inputs, intermediate expression/cell matrices, response maps and private protocol records are not shipped. The filename private_dependencies/unavailable_protocol.txt is an explicit unavailable-input marker in sanitized historical code, not a request to copy private instructions into a public checkout. Do not fabricate that file or disable hash guards. A reviewed public replacement and corresponding non-scientific provenance adaptation are still required before those initializer paths can execute.

## Known execution gaps

Some initializers import private authorization artifacts; source locks also bind historical paths, scripts or unshipped protocol records. Sanitizing paths changes script hashes, so a historical self-hash cannot be assumed to match the public copy. Original and sanitized code hashes are both recorded. Some mixed clinical source files retain excluded exploratory-cohort branches. Active numerical producers were supplemented where the supplied map pointed only to drawing/assembly, but comprehensive producer-to-output lineage and external-input closure are not yet demonstrated.

Public deployment is NOT_READY until these documented portability and input-contract gaps are resolved without changing calculations, scientific locks or frozen results. This is stronger than a missing raw-data download: some required historical records cannot appropriately be redistributed. The candidate still provides inspectable code, exact definitions, original derived summaries and explicit limits.

## Environments

Historical version records are module-specific. No universal environment lockfile is claimed. The historical SciPy version for the topology module was not recoverable. Installing a current version is not evidence of the historical version or proof of numerical equivalence.

## Manifest boundary

release/RELEASE_MANIFEST.tsv covers payload files; it excludes itself and release/SHA256SUMS.txt to avoid circular hashes. SHA256SUMS.txt includes every repository payload file and RELEASE_MANIFEST.tsv, excluding only itself. Git internals and ignored local workspaces are not release payload. Source hashes describe the original artifact; public hashes describe the sanitized artifact. No remote, push or public release is created by these instructions.
