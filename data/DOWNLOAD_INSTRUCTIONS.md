# Obtaining source data

1. Find the accession, original study and intended role in SOURCE_DATA_MANIFEST.tsv.
2. For a GEO series, open its linked series page, review the study publication and series/sample metadata, and obtain the relevant supplementary matrix, annotation and processed files. Do not substitute a different accession or infer missing response maps.
3. For ERP117672, follow the linked ENA study and Hong et al. publication for the molecular files and clinical mapping. The 40-sample molecular population and 35-patient response-evaluable population are different.
4. For Wu ST-CODEX, use the exact Zenodo record identified in the manifest; obtain the released spatial expression, annotations, coordinates and registered CODEX inputs needed for the chosen module. Adjacent sections are not the same cells.
5. Review original terms and retain source provenance/checksums locally. Do not upload downloaded files with this repository. No large download is required for package verification.
6. Follow the file references and locked input contracts for the chosen module. File-by-file acquisition and preparation are not fully resolved for every module; read docs/REPRODUCIBILITY_NOTES.md before attempting an analysis. Do not bypass source hashes or create substitute clinical mappings.

Source download instructions are accession-level guidance, not a verified automatic acquisition pipeline. Original studies and source terms remain authoritative.
