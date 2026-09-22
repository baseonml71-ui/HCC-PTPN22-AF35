from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


PROJECT = Path(r"./HCC_ICI_project")
META = PROJECT / "01_metadata" / "phase19b_published_benchmark"
TABLES = PROJECT / "05_results" / "tables" / "phase19b_published_benchmark"
PROCESSED = PROJECT / "03_processed_data" / "phase19b_published_benchmark"

PATHS = {
    "candidate_pool": META / "PHASE19B_PUBLISHED_SIGNATURE_CANDIDATE_POOL.tsv",
    "contract": META / "PHASE19B_EXACT_PUBLISHED_SIGNATURE_CONTRACT.tsv",
    "tier": META / "PHASE19B_PUBLISHED_SIGNATURE_TIER_ASSIGNMENT.tsv",
    "erp_scores": TABLES / "ERP117672_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv",
    "g302_scores": TABLES / "GSE302495_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv",
    "relationships": TABLES / "AF35_PUBLISHED_SIGNATURE_MOLECULAR_RELATIONSHIPS.tsv",
    "blind_object": PROCESSED / "PHASE19B_RESPONSE_BLIND_OBJECT.rds",
    "blind_freeze": META / "PHASE19B_RESPONSE_BLIND_FREEZE.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


for path in PATHS.values():
    if not path.is_file():
        raise FileNotFoundError(path)

candidate = read_tsv(PATHS["candidate_pool"])
contract = read_tsv(PATHS["contract"])
tier = read_tsv(PATHS["tier"])
erp = read_tsv(PATHS["erp_scores"])
g302 = read_tsv(PATHS["g302_scores"])
relationships = read_tsv(PATHS["relationships"])
blind_freeze = json.loads(PATHS["blind_freeze"].read_text(encoding="utf-8"))

expected_ids = ["S005", "S006", "S007", "S008", "S009", "S010", "S011", "S012", "S013", "S014", "S015", "S021", "S022", "S024", "P18R13623"]
if [row["paper_id"] for row in candidate] != expected_ids:
    raise RuntimeError("Candidate pool does not match the frozen Phase18 order")
if [row["paper_id"] for row in contract] != expected_ids:
    raise RuntimeError("Contract does not match the frozen candidate pool")
eligible = [row for row in contract if row["exact_definition_status"] == "ELIGIBLE_TIER_A"]
no_go = [row for row in contract if row["exact_definition_status"] == "PUBLISHED_SIGNATURE_DEFINITION_NO_GO"]
if len(eligible) != 1 or eligible[0]["paper_id"] != "S007" or len(no_go) != 14:
    raise RuntimeError("Eligible/no-go set is not frozen as expected")
if len(eligible[0]["gene_list"].split(";")) != 20:
    raise RuntimeError("S007 exact gene list is not 20 genes")
if len(erp) != 40 or len(g302) != 38:
    raise RuntimeError("Response-blind score populations changed")
for cohort_rows in (erp, g302):
    if any(row["outcome_joined"].upper() != "FALSE" for row in cohort_rows):
        raise RuntimeError("Outcome label joined before lock")
    forbidden = {"response", "response_group", "group", "best_response", "cutpoint", "roc", "auroc"}
    if forbidden.intersection(cohort_rows[0]):
        raise RuntimeError("Forbidden outcome or cutoff column in response-blind score table")
if len(relationships) != 1 or relationships[0]["outcome_joined"].upper() != "FALSE":
    raise RuntimeError("Molecular relationship table is not response-blind")
if relationships[0]["near_duplicate_predefined_flag"].upper() != "FALSE":
    raise RuntimeError("Unexpected near-duplicate flag")
if [row["paper_id"] for row in tier if row["tier"] == "A"] != ["S007"]:
    raise RuntimeError("Tier A set is not frozen as expected")
if blind_freeze["blind_object_sha256"] != sha256(PATHS["blind_object"]):
    raise RuntimeError("Response-blind RDS hash mismatch")

lock = {
    "created_at": datetime.now().astimezone().isoformat(),
    "status": "PHASE19B_RESPONSE_BLIND_TABLES_LOCKED_BEFORE_OUTCOME_JOIN",
    "outcome_joined": False,
    "response_map_read_by_scoring_stage": False,
    "eligible_signature_n": 1,
    "eligible_signature": "S007_ISG_20",
    "eligible_peer_reviewed_tier_a_paper_n": 1,
    "no_go_row_n": 14,
    "gate_prerequisite": "FAIL_LT_2_INDEPENDENT_PEER_REVIEWED_TIER_A",
    "file_sha256": {name: sha256(path) for name, path in PATHS.items()},
}
lock_path = META / "PHASE19B_RESPONSE_BLIND_TABLE_LOCK.json"
lock_path.write_text(json.dumps(lock, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"status": "PASS", "lock": str(lock_path), "lock_sha256": sha256(lock_path)}, ensure_ascii=False))
