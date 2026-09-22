import csv
import gzip
import hashlib
import json
import math
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "10_software" / "phase3_py"))

import numpy as np


PROJECTION_FILE = PROJECT_ROOT / "03_processed_data/phase19a_af35_clonotype/GSE235863_AF35_RESPONSE_BLIND_PROJECTION.npz"
FREEZE_FILE = PROJECT_ROOT / "01_metadata/phase19a_af35_clonotype/PHASE19A_AF35_RESPONSE_BLIND_PROJECTION_FREEZE.json"
MAPPING_FILE = PROJECT_ROOT / "03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz"
FATE_FILE = PROJECT_ROOT / "05_results/tables/temporal/GSE235863_BASELINE_CLONE_FATE_TABLE.tsv"
WORK_DIR = PROJECT_ROOT / "10_intermediate_files/phase19"
TABLE_JSON = WORK_DIR / "phase19a_table_data.json"
SUMMARY_JSON = WORK_DIR / "phase19a_summary.json"
LOG_FILE = PROJECT_ROOT / "09_logs/phase19/02_phase19a_patient_analysis.log"
SEED = 20260831


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def truth(value):
    return str(value).upper() == "TRUE"


def rank_average(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(values.size, dtype=float)
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.unique(x[ok]).size < 2 or np.unique(y[ok]).size < 2:
        return math.nan
    return float(np.corrcoef(rank_average(x[ok]), rank_average(y[ok]))[0, 1])


def summarize_effects(effects, rng, draws=10000):
    values = np.asarray([value for value in effects.values() if np.isfinite(value)], dtype=float)
    patients = [patient for patient, value in effects.items() if np.isfinite(value)]
    if values.size == 0:
        return {
            "eligible_patient_n": 0, "mean_patient_effect": None, "median_patient_effect": None,
            "bootstrap_ci_low": None, "bootstrap_ci_high": None, "positive_patient_n": 0,
            "negative_patient_n": 0, "zero_patient_n": 0, "positive_direction_fraction": None,
            "loo_mean_min": None, "loo_mean_max": None, "loo_positive_direction_fraction": None,
        }
    boot = np.array([rng.choice(values, size=values.size, replace=True).mean() for _ in range(draws)])
    loo = np.array([np.delete(values, i).mean() for i in range(values.size)]) if values.size > 1 else np.array([values[0]])
    return {
        "eligible_patient_n": int(values.size),
        "eligible_patients": ";".join(patients),
        "mean_patient_effect": float(values.mean()),
        "median_patient_effect": float(np.median(values)),
        "bootstrap_ci_low": float(np.quantile(boot, 0.025)),
        "bootstrap_ci_high": float(np.quantile(boot, 0.975)),
        "bootstrap_draws": draws,
        "positive_patient_n": int(np.sum(values > 0)),
        "negative_patient_n": int(np.sum(values < 0)),
        "zero_patient_n": int(np.sum(values == 0)),
        "positive_direction_fraction": float(np.mean(values > 0)),
        "loo_mean_min": float(loo.min()),
        "loo_mean_max": float(loo.max()),
        "loo_positive_direction_fraction": float(np.mean(loo > 0)),
    }


def group_effect(rows, flag_key, score_key):
    one = [row[score_key] for row in rows if row[flag_key] is True]
    zero = [row[score_key] for row in rows if row[flag_key] is False]
    if not one or not zero:
        return math.nan, len(one), len(zero), math.nan, math.nan
    mean_one = float(np.mean(one))
    mean_zero = float(np.mean(zero))
    return mean_one - mean_zero, len(one), len(zero), mean_one, mean_zero


def state_stratified_effect(rows, flag_key, score_key, state_key="sub_cluster"):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row[state_key]].append(row)
    effects = []
    weights = []
    for state_rows in grouped.values():
        effect, n1, n0, _, _ = group_effect(state_rows, flag_key, score_key)
        if np.isfinite(effect):
            effects.append(effect)
            weights.append(min(n1, n0))
    if not weights:
        return math.nan
    return float(np.average(np.asarray(effects), weights=np.asarray(weights)))


def categorical_stratified_effect(rows, flag_key, score_key, category_key):
    return state_stratified_effect(rows, flag_key, score_key, category_key)


def equal_cell_downsample_effect(rows, flag_key, score_key, rng, draws=500):
    one = np.asarray([row[score_key] for row in rows if row[flag_key] is True], dtype=float)
    zero = np.asarray([row[score_key] for row in rows if row[flag_key] is False], dtype=float)
    take = min(one.size, zero.size)
    if take == 0:
        return math.nan
    values = [rng.choice(one, take, replace=False).mean() - rng.choice(zero, take, replace=False).mean() for _ in range(draws)]
    return float(np.mean(values))


def patient_rows_with_aggregate(analysis, context_one, context_zero, rows_by_patient, flag_key, score_key, rng):
    output = []
    effects = {}
    for patient in sorted(rows_by_patient):
        rows = rows_by_patient[patient]
        effect, n1, n0, mean1, mean0 = group_effect(rows, flag_key, score_key)
        if not np.isfinite(effect):
            continue
        response = rows[0]["response_group"]
        effects[patient] = effect
        output.append({
            "analysis": analysis,
            "row_type": "PATIENT",
            "patient_id": patient,
            "response_group": response,
            "context_positive": context_one,
            "context_reference": context_zero,
            "n_positive_cells": n1,
            "n_reference_cells": n0,
            "mean_positive_af35": mean1,
            "mean_reference_af35": mean0,
            "patient_effect_positive_minus_reference": effect,
            "score_projection": score_key,
            "biological_replicate": "patient",
        })
    aggregate = summarize_effects(effects, rng)
    output.append({
        "analysis": analysis,
        "row_type": "AGGREGATE",
        "patient_id": "AGGREGATE",
        "response_group": "ALL",
        "context_positive": context_one,
        "context_reference": context_zero,
        "score_projection": score_key,
        "biological_replicate": "patient",
        **aggregate,
    })
    return output, effects, aggregate


def sensitivity_row(module, analysis, method, effects, rng, note):
    return {
        "module": module,
        "analysis": analysis,
        "sensitivity_method": method,
        **summarize_effects(effects, rng),
        "aggregate_direction": "POSITIVE" if np.nanmean(list(effects.values())) > 0 else "NONPOSITIVE",
        "note": note,
        "biological_replicate": "patient",
    }


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    messages = []

    def record(message):
        text = str(message)
        messages.append(text)
        print(text, flush=True)

    record("Phase 19A patient-aware analysis")
    record(f"Started: {datetime.now().astimezone().isoformat()}")
    freeze = json.loads(FREEZE_FILE.read_text(encoding="utf-8"))
    if freeze["outcome_joined"] is not False or freeze["usable_gene_count"] < 28:
        raise RuntimeError("Response-blind projection freeze is not eligible")
    if sha256(PROJECTION_FILE) != freeze["npz_sha256"]:
        raise RuntimeError("Projection NPZ hash mismatch")

    projection = np.load(PROJECTION_FILE, allow_pickle=False)
    score_map = {
        barcode: (float(primary), float(rank))
        for barcode, primary, rank in zip(
            projection["barcode"].tolist(), projection["af35_primary"], projection["af35_rank_sensitivity"]
        )
    }

    rows = []
    with gzip.open(MAPPING_FILE, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for raw in reader:
            if raw["barcode"] not in score_map:
                raise RuntimeError(f"Projection missing barcode {raw['barcode']}")
            primary, rank = score_map[raw["barcode"]]
            rows.append({
                "patient": raw["patient"],
                "response_group": raw["response_group"],
                "sample": raw["sample"],
                "timepoint": raw["timepoint"],
                "tissue": raw["tissue"],
                "barcode": raw["barcode"],
                "clone_id": raw["clone.id"],
                "clonotype_key": raw["clonotype_key"],
                "cdr3_pair_aa": raw["cdr3_pair_aa"],
                "clone_size": int(raw["clone.size"]),
                "major_cluster": raw["major_cluster"],
                "sub_cluster": raw["sub_cluster"],
                "persistence_eligible": truth(raw["persistence_eligible"]),
                "persistent": truth(raw["persistent_within_tissue"]),
                "expanded": int(raw["clone.size"]) >= 2,
                "af35_primary": primary,
                "af35_rank": rank,
            })
    if len(rows) != 58872:
        raise RuntimeError("Mapped row count changed")

    fate = {}
    fate_rows = []
    with FATE_FILE.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for raw in reader:
            if raw["analysis_role"] != "PRIMARY_BLOOD":
                continue
            key = (raw["patient"], raw["clonotype_key"])
            item = {
                "future_observed": truth(raw["future_persistent_observed"]),
                "baseline_size_category": raw["baseline_clone_size_category"],
                "future_abundance_change": as_float(raw["future_abundance_change_persistent"]),
            }
            fate[key] = item
            fate_rows.append({**raw, **item})

    for row in rows:
        item = fate.get((row["patient"], row["clonotype_key"]))
        row["future_observed"] = item["future_observed"] if item else None
        row["baseline_size_category"] = item["baseline_size_category"] if item else ""

    expanded_by_patient = defaultdict(list)
    persistent_by_patient = defaultdict(list)
    temporal_by_patient = defaultdict(list)
    for row in rows:
        expanded_by_patient[row["patient"]].append(row)
        if row["persistence_eligible"]:
            persistent_by_patient[row["patient"]].append(row)
        if row["timepoint"] == "pre" and row["tissue"] == "P" and row["future_observed"] is not None:
            temporal_by_patient[row["patient"]].append(row)

    expanded_table, expanded_effects, expanded_summary = patient_rows_with_aggregate(
        "EXPANDED_VS_SINGLETON", "expanded", "singleton", expanded_by_patient, "expanded", "af35_primary", rng
    )
    persistent_table, persistent_effects, persistent_summary = patient_rows_with_aggregate(
        "PERSISTENT_VS_NONPERSISTENT", "persistent", "nonpersistent", persistent_by_patient, "persistent", "af35_primary", rng
    )
    temporal_table, temporal_effects, temporal_summary = patient_rows_with_aggregate(
        "BASELINE_FUTURE_OBSERVED_VS_NOT_OBSERVED", "future_observed", "not_observed", temporal_by_patient,
        "future_observed", "af35_primary", rng
    )

    sensitivities = []
    for analysis, by_patient, flag, primary_effects in [
        ("EXPANDED_VS_SINGLETON", expanded_by_patient, "expanded", expanded_effects),
        ("PERSISTENT_VS_NONPERSISTENT", persistent_by_patient, "persistent", persistent_effects),
    ]:
        sensitivities.append(sensitivity_row("CLONOTYPE_CONTEXT", analysis, "PRIMARY_GENEWISE_Z", primary_effects, rng, "Primary frozen AF35 projection"))
        rank_effects = {p: group_effect(rs, flag, "af35_rank")[0] for p, rs in by_patient.items()}
        sensitivities.append(sensitivity_row("CLONOTYPE_CONTEXT", analysis, "RANK_STANDARDIZED_PROJECTION", rank_effects, rng, "Response-blind rank sensitivity"))
        for major in ("CD4", "CD8"):
            effects = {p: group_effect([r for r in rs if r["major_cluster"] == major], flag, "af35_primary")[0] for p, rs in by_patient.items()}
            effects = {p: value for p, value in effects.items() if np.isfinite(value)}
            sensitivities.append(sensitivity_row("CLONOTYPE_CONTEXT", analysis, f"{major}_ONLY", effects, rng, "Frozen author major-cluster label"))
        state_effects = {p: state_stratified_effect(rs, flag, "af35_primary") for p, rs in by_patient.items()}
        state_effects = {p: value for p, value in state_effects.items() if np.isfinite(value)}
        sensitivities.append(sensitivity_row("CLONOTYPE_CONTEXT", analysis, "AUTHOR_STATE_STRATIFIED", state_effects, rng, "Exact author subcluster; min-group-cell weighted"))
        down_effects = {p: equal_cell_downsample_effect(rs, flag, "af35_primary", rng) for p, rs in by_patient.items()}
        down_effects = {p: value for p, value in down_effects.items() if np.isfinite(value)}
        sensitivities.append(sensitivity_row("CLONOTYPE_CONTEXT", analysis, "EQUAL_CELL_NUMBER_DOWNSAMPLING", down_effects, rng, "500 fixed-seed draws within patient/context"))

    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "PRIMARY_GENEWISE_Z", temporal_effects, rng, "Exact frozen primary blood set"))
    rank_effects = {p: group_effect(rs, "future_observed", "af35_rank")[0] for p, rs in temporal_by_patient.items()}
    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "RANK_STANDARDIZED_PROJECTION", rank_effects, rng, "Response-blind rank sensitivity"))
    size_effects = {p: categorical_stratified_effect(rs, "future_observed", "af35_primary", "baseline_size_category") for p, rs in temporal_by_patient.items()}
    size_effects = {p: value for p, value in size_effects.items() if np.isfinite(value)}
    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "BASELINE_CLONE_SIZE_STRATIFIED", size_effects, rng, "Frozen 1, 2-4, >=5 baseline clone-size categories"))
    state_effects = {p: state_stratified_effect(rs, "future_observed", "af35_primary") for p, rs in temporal_by_patient.items()}
    state_effects = {p: value for p, value in state_effects.items() if np.isfinite(value)}
    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "AUTHOR_STATE_STRATIFIED", state_effects, rng, "Exact author subcluster; min-group-cell weighted"))

    followup_by_patient = defaultdict(list)
    for row in rows:
        if row["timepoint"] == "post" and row["tissue"] == "P" and row["patient"] in temporal_by_patient:
            followup_by_patient[row["patient"]].append(row)
    depth_effects = {}
    for patient, baseline_rows in temporal_by_patient.items():
        followup = followup_by_patient[patient]
        if len(followup) < 580:
            continue
        draw_effects = []
        for _ in range(500):
            sampled = rng.choice(len(followup), size=580, replace=False)
            observed = {followup[i]["clonotype_key"] for i in sampled}
            one = [row["af35_primary"] for row in baseline_rows if row["clonotype_key"] in observed]
            zero = [row["af35_primary"] for row in baseline_rows if row["clonotype_key"] not in observed]
            if one and zero:
                draw_effects.append(float(np.mean(one) - np.mean(zero)))
        if draw_effects:
            depth_effects[patient] = float(np.median(draw_effects))
    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "FOLLOWUP_COMMON_DEPTH_580", depth_effects, rng, "500 fixed-seed technical draws; follow-up absence remains non-detection"))

    cdr3_effects = {}
    followup_cdr3 = {p: {row["cdr3_pair_aa"] for row in rs} for p, rs in followup_by_patient.items()}
    for patient, baseline_rows in temporal_by_patient.items():
        copied = []
        for row in baseline_rows:
            item = dict(row)
            item["future_observed_cdr3"] = row["cdr3_pair_aa"] in followup_cdr3[patient]
            copied.append(item)
        cdr3_effects[patient] = group_effect(copied, "future_observed_cdr3", "af35_primary")[0]
    sensitivities.append(sensitivity_row("TEMPORAL", "BASELINE_FUTURE_OBSERVED", "PAIRED_CDR3_ALPHA_BETA_IDENTITY", cdr3_effects, rng, "Patient + paired CDR3alpha + paired CDR3beta"))

    response_effects = defaultdict(list)
    for patient, effect in temporal_effects.items():
        response_effects[temporal_by_patient[patient][0]["response_group"]].append(effect)
    for response_group, values in sorted(response_effects.items()):
        sensitivities.append({
            "module": "TEMPORAL", "analysis": "BASELINE_FUTURE_OBSERVED", "sensitivity_method": f"RESPONSE_{response_group}_DESCRIPTIVE",
            "eligible_patient_n": len(values), "mean_patient_effect": float(np.mean(values)), "median_patient_effect": float(np.median(values)),
            "positive_patient_n": int(np.sum(np.asarray(values) > 0)), "positive_direction_fraction": float(np.mean(np.asarray(values) > 0)),
            "aggregate_direction": "POSITIVE" if np.mean(values) > 0 else "NONPOSITIVE", "note": "Descriptive only; no interaction test", "biological_replicate": "patient",
        })

    baseline_clone_scores = defaultdict(list)
    for patient, baseline_rows in temporal_by_patient.items():
        clone_values = defaultdict(list)
        for row in baseline_rows:
            clone_values[row["clonotype_key"]].append(row["af35_primary"])
        for clone_key, values in clone_values.items():
            change = fate[(patient, clone_key)]["future_abundance_change"]
            if np.isfinite(change):
                baseline_clone_scores[patient].append((float(np.mean(values)), change))
    abundance_rho = {p: spearman([x for x, _ in pairs], [y for _, y in pairs]) for p, pairs in baseline_clone_scores.items()}
    abundance_rho = {p: value for p, value in abundance_rho.items() if np.isfinite(value)}
    sensitivities.append(sensitivity_row("TEMPORAL", "FUTURE_ABUNDANCE_CHANGE_PERSISTENT_ONLY", "FROZEN_V6_OUTCOME_TABLE_SPEARMAN", abundance_rho, rng, "Cannot rescue primary detection result"))

    tables = {
        "PHASE19A_AF35_GENE_PROJECTION_PARAMETERS.tsv": freeze["genes"],
        "AF35_EXPANDED_SINGLETON_PATIENT_EFFECTS.tsv": expanded_table,
        "AF35_PERSISTENT_NONPERSISTENT_PATIENT_EFFECTS.tsv": persistent_table,
        "AF35_BASELINE_TO_FUTURE_OBSERVED_PATIENT_EFFECTS.tsv": temporal_table,
        "AF35_CLONOTYPE_TEMPORAL_SENSITIVITY.tsv": sensitivities,
    }
    TABLE_JSON.write_text(json.dumps(tables, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    major_sensitivity_rows = [row for row in sensitivities if row["sensitivity_method"] in {
        "RANK_STANDARDIZED_PROJECTION", "AUTHOR_STATE_STRATIFIED", "EQUAL_CELL_NUMBER_DOWNSAMPLING",
        "BASELINE_CLONE_SIZE_STRATIFIED", "FOLLOWUP_COMMON_DEPTH_580", "PAIRED_CDR3_ALPHA_BETA_IDENTITY",
    }]
    all_major_positive = all(row["aggregate_direction"] == "POSITIVE" for row in major_sensitivity_rows)
    identity_row = next(row for row in sensitivities if row["sensitivity_method"] == "PAIRED_CDR3_ALPHA_BETA_IDENTITY")
    if (
        expanded_summary["positive_direction_fraction"] >= 0.8
        and persistent_summary["positive_direction_fraction"] >= 0.8
        and temporal_summary["positive_direction_fraction"] >= 0.8
        and all_major_positive
        and identity_row["aggregate_direction"] == "POSITIVE"
    ):
        gate = "AF35_CLONOTYPE_TEMPORAL_BRIDGE_SUPPORTED"
    elif any(item["positive_direction_fraction"] >= 0.8 for item in (expanded_summary, persistent_summary, temporal_summary)):
        gate = "AF35_CLONOTYPE_TEMPORAL_BRIDGE_PARTIAL"
    else:
        gate = "AF35_CLONOTYPE_TEMPORAL_BRIDGE_NOT_SUPPORTED"

    summary = {
        "created_at": datetime.now().astimezone().isoformat(),
        "seed": SEED,
        "projection_freeze_sha256": sha256(FREEZE_FILE),
        "projection_npz_sha256": sha256(PROJECTION_FILE),
        "mapped_cell_n": len(rows),
        "expanded": expanded_summary,
        "persistent": persistent_summary,
        "temporal": temporal_summary,
        "major_sensitivities_positive": all_major_positive,
        "clone_identity_temporal_positive_fraction": identity_row["positive_direction_fraction"],
        "gate": gate,
        "cellular_source_audit": "NOT_RUN_RESOURCE_NOT_FROZEN",
        "ptpn22_af35_pattern": None,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    record(f"Expanded positive: {expanded_summary['positive_patient_n']}/{expanded_summary['eligible_patient_n']}")
    record(f"Persistent positive: {persistent_summary['positive_patient_n']}/{persistent_summary['eligible_patient_n']}")
    record(f"Temporal positive: {temporal_summary['positive_patient_n']}/{temporal_summary['eligible_patient_n']}")
    record(f"Gate: {gate}")
    record(f"Completed: {datetime.now().astimezone().isoformat()}")
    LOG_FILE.write_text("\n".join(messages) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
