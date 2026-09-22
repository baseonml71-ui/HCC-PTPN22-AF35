# Public repository release gate

**Gate: NOT_READY**

Local review status: HCC_PUBLIC_REPOSITORY_V1_RC1_AUTHOR_REVIEW.

The candidate is sanitized for local author review, but unresolved public execution contracts prevent a publication-ready claim. It has not been pushed, uploaded or published. Scientific results are unchanged.

| Requirement | Assessment |
|---|---|
| AF35 definition integrity | PASS: 35 unique ordered genes match recovered vector; equal-weight, unsigned, PTPN22 excluded |
| Active-script completeness | NOT COMPLETE: supplied labels corrected; numerical producers/helpers included; full dynamic input and producer-output closure remains unresolved |
| Local-path sanitization | PASS for release payload: project roots replaced by relative roots; private local resources replaced with explicit unavailable-input markers |
| Credential/privacy scan | No identified unresolved credentials, private tokens, direct PHI or unintended absolute local paths in admitted payload; details in security audit |
| Third-party asset safety | PASS by exclusion: no figure binaries or decorative liver illustration; every admitted non-code asset inventoried |
| Raw-data redistribution | PASS: no external raw or individual-level data included; accessions and source terms supplied |
| Environment documentation | PASS as historical records; no complete lockfile; historical topology SciPy version explicitly unavailable |
| Citation metadata | PASS: official CFF 1.2.0 schema; six authors; two supplied ORCIDs; publication URL/DOI pending |
| MIT license | PASS: standard text with author-specified copyright; external materials retain their terms |
| README correctness | Static package checks pass; reproduction limitations and deliberately failing prerequisite check are explicit |
| Release manifest | Payload inventory excludes itself and checksums to avoid a circular hash; checksums include manifest |
| Checksums | SHA256SUMS covers all payload except itself; verify_release.py checks both set equality and content |

## Actions before public publication

1. Resolve the private protocol-input dependencies listed in provenance/EXECUTION_GAPS.tsv with a reviewed public substitute or non-scientific adapter. Do not publish private instructions, fabricate historical records or bypass frozen scientific checks.
2. Complete module-level acquisition/preparation contracts and current producer-to-output mapping, including dynamic helper inputs and historical self-hash/path dependencies. Dataset accessions alone do not resolve these dependencies.
3. Separate the current clinical path from historical branches in shared source entrypoints, preserving calculations and recording an auditable correspondence.
4. Repeat the bounded static/integrity/privacy checks after those changes. Do not claim a successful raw-data rerun unless one is separately authorized and actually completed.
5. After the gate passes, obtain author authorization to create GitHub hosting and Zenodo archival. Fill the actual repository URL and issued DOI then, not in advance.

The missing historical SciPy version remains a documented limitation, not a fabricated pin. Figure 3A rights do not block this asset-free candidate because the illustration and composites are absent. No newly identified credential or PHI issue remains; the blocking issues concern execution contracts and scientific-code lineage.
