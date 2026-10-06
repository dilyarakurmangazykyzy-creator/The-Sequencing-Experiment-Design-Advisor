#!/usr/bin/env python3
"""Deterministic, transparent sequencing-design planning model (Python standard library).

Use: python advisor.py --input config/scenario.json
     python advisor.py --organism bacteria --application variant --budget 1200 --samples 6
The empirical pilot measures anchor-supported reference breadth. It is never relabelled
as variant recall, assembly completeness, RNA-seq power or clinical sensitivity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
MODEL_VERSION = "1.0.0"
DEFAULTS: dict[str, Any] = {
    "organism": "bacteria", "application": "variant", "budget_usd": 1200,
    "sample_count": 6, "allele_fraction": 1, "detection_probability": 0.95,
    "min_alt_reads": 3, "mapping_fraction": 0.95, "usable_fraction": 0.85,
    "repeat_span_kb": 8, "replicates_per_condition": 3, "condition_count": 2,
    "expression_pairs_millions": 30, "expression_long_reads_millions": 1,
    "procurement": "shared", "platform": "auto", "assembly_goal": "resolve_repeats",
    "expression_goal": "gene_counts",
}
GENOME_MB = {"bacteria": 4.6, "yeast": 12.1, "human": 3100, "custom": 100}
ENUMS = {
    "organism": set(GENOME_MB), "application": {"variant", "expression", "assembly"},
    "procurement": {"shared", "dedicated"},
    "platform": {"auto", "illumina", "pacbio", "nanopore"},
    "assembly_goal": {"draft", "resolve_repeats"},
    "expression_goal": {"gene_counts", "full_length_isoforms"},
}
NUMBERS = {
    "budget_usd": (1, 1e9, False), "sample_count": (1, 100000, True),
    "genome_size_mb": (0.001, 1e6, False), "allele_fraction": (0.00001, 1, False),
    "detection_probability": (0.5, 0.99999, False), "min_alt_reads": (1, 100, True),
    "mapping_fraction": (0.01, 1, False), "usable_fraction": (0.01, 1, False),
    "repeat_span_kb": (0, 10000, False), "replicates_per_condition": (1, 10000, True),
    "condition_count": (1, 1000, True), "expression_pairs_millions": (0.01, 10000, False),
    "expression_long_reads_millions": (0.01, 1000, False),
}


class InputError(ValueError):
    """A user-fixable scenario or configuration error."""


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    selected = Path(path) if path is not None else ROOT / "config" / "costs.json"
    try:
        config = read_json(selected)
    except (OSError, ValueError) as error:
        raise InputError(f"Cannot read cost configuration: {error}") from error
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise InputError("Cost configuration must be a schema_version=1 JSON object.")
    required_global = ("uncertainty_fraction", "contingency_fraction", "compute_hour_usd",
                       "storage_gb_month_usd", "storage_months", "qc_per_sample_usd",
                       "analysis_setup_hours", "expression_min_pairs_millions",
                       "expression_min_long_reads_millions", "expression_count_threshold")
    required_platform = ("run_yield_gb", "run_hours", "reagent_run_usd",
                         "shared_sequencing_usd_per_gb", "shared_min_sequencing_usd",
                         "library_per_sample_usd", "rna_library_per_sample_usd",
                         "instrument_usd_per_hour", "compute_hours_per_gb",
                         "assembly_compute_multiplier", "stored_gb_per_sequence_gb", "read_length_bp")
    for key in required_global:
        _config_number(config, key, allow_zero=True)
    if config["uncertainty_fraction"] >= 1:
        raise InputError("uncertainty_fraction must be less than 1.")
    for platform in ("illumina", "pacbio", "nanopore"):
        try:
            details = config["platforms"][platform]
        except (KeyError, TypeError) as error:
            raise InputError(f"Missing platform configuration: {platform}") from error
        for key in required_platform:
            _config_number(details, key, allow_zero=key not in (
                "run_yield_gb", "run_hours", "shared_sequencing_usd_per_gb", "read_length_bp"))
        for floor_group in ("assembly_floor_x",):
            if not isinstance(config.get(floor_group), dict):
                raise InputError(f"Missing {floor_group} configuration.")
            _config_number(config[floor_group], platform)
    for organism in GENOME_MB:
        if not isinstance(config.get("variant_floor_x"), dict):
            raise InputError("Missing variant_floor_x configuration.")
        _config_number(config["variant_floor_x"], organism)
    threshold = config["expression_count_threshold"]
    if int(threshold) != threshold or not 1 <= threshold <= 100:
        raise InputError("expression_count_threshold must be an integer from 1 to 100.")
    return config


def _config_number(obj: dict, key: str, allow_zero: bool = False) -> None:
    value = obj.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise InputError(f"Configuration '{key}' must be a finite number.")
    if value < 0 or (value == 0 and not allow_zero):
        raise InputError(f"Configuration '{key}' must be {'nonnegative' if allow_zero else 'positive'}.")


def normalize_input(scenario: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(scenario, dict):
        raise InputError("Scenario must be a JSON object.")
    unknown = set(scenario) - set(DEFAULTS) - {"genome_size_mb"}
    if unknown:
        raise InputError("Unknown scenario field(s): " + ", ".join(sorted(unknown)))
    normalized = {**DEFAULTS, **scenario}
    for key, values in ENUMS.items():
        if not isinstance(normalized[key], str) or normalized[key] not in values:
            raise InputError(f"{key} must be one of: {', '.join(sorted(values))}.")
    if "genome_size_mb" not in normalized:
        normalized["genome_size_mb"] = GENOME_MB[normalized["organism"]]
    for key, (lower, upper, integer) in NUMBERS.items():
        value = normalized[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise InputError(f"{key} must be a finite number.")
        if not lower <= value <= upper:
            raise InputError(f"{key} must be between {lower:g} and {upper:g}.")
        if integer and int(value) != value:
            raise InputError(f"{key} must be a whole number.")
        normalized[key] = int(value) if integer else float(value)
    return normalized


def poisson_tail(mean: float, min_reads: int) -> float:
    """P(X >= min_reads) for X~Poisson(mean); stable in the supported range."""
    if mean <= 0:
        return 0.0
    if min_reads < 1:
        return 1.0
    logs = [-mean + i * math.log(mean) - math.lgamma(i + 1) for i in range(min_reads)]
    maximum = max(logs)
    cdf = math.exp(maximum) * math.fsum(math.exp(value - maximum) for value in logs)
    return min(1.0, max(0.0, 1.0 - cdf))


def required_poisson_mean(probability: float, min_reads: int) -> float:
    lower, upper = 0.0, float(max(1, min_reads))
    while poisson_tail(upper, min_reads) < probability:
        upper *= 2
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if poisson_tail(midpoint, min_reads) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return upper


def load_evidence(root: Path | str | None = None) -> dict[str, Any]:
    base = Path(root) if root is not None else ROOT
    evidence: dict[str, Any] = {
        "status": "absent", "scope": "Limited bacterial reference-mapping pilot",
        "endpoint": "anchor_supported_reference_breadth", "saturation": None, "summary": None,
        "issues": [], "recommendation_use": "Diagnostic comparison only; cannot validate calling depth.",
        "limitations": [
            "Anchor-supported mapped breadth is not variant recall or an established mapping-quality score.",
            "A plateau below complete breadth may reflect reference mismatch, repeats or aligner limits.",
            "A bacterial pilot cannot establish human, RNA-seq or assembly depth requirements.",
            "Repeated downsampling seeds measure sampling variability, not biological variability.",
        ],
    }
    for name in ("saturation", "summary"):
        path = base / "results" / f"{name}.json"
        if not path.exists():
            evidence["issues"].append(f"results/{name}.json is not available; run the analysis workflow.")
            continue
        try:
            item = read_json(path)
            if not isinstance(item, dict):
                raise ValueError("expected a JSON object")
            evidence[name] = item
        except (OSError, ValueError) as error:
            evidence["issues"].append(f"Cannot load {name} evidence: {error}")
    saturation = evidence["saturation"]
    if saturation is not None:
        supported_endpoints = ("anchor_supported_reference_breadth", "reference_breadth_primary_alignments_mapq_ge_20")
        if saturation.get("endpoint") not in supported_endpoints:
            evidence["issues"].append("Unrecognized saturation endpoint; recommendation does not use its values.")
        else:
            platforms = saturation.get("platforms", {})
            valid = isinstance(platforms, dict) and any(
                isinstance(dataset, dict) and bool(dataset.get("points"))
                for dataset in platforms.values())
            if valid:
                evidence["status"] = "limited_pilot"
                evidence["endpoint"] = saturation["endpoint"]
                if saturation["endpoint"] == "reference_breadth_primary_alignments_mapq_ge_20":
                    evidence["limitations"][0] = "Reference breadth of primary alignments with MAPQ ≥20 is not variant recall, expression power or assembly completeness."
            else:
                evidence["issues"].append("No measured saturation points are present.")
    if evidence["status"] == "absent" and evidence["summary"] is not None:
        evidence["status"] = "incomplete"
    return evidence


def cost_for_raw_gb(raw_gb: float, scenario: dict, platform: str, config: dict) -> dict:
    p = config["platforms"][platform]
    n = scenario["sample_count"]
    application = scenario["application"]
    run_equivalents = raw_gb / p["run_yield_gb"]
    dedicated = scenario["procurement"] == "dedicated"
    billed_runs = math.ceil(run_equivalents - 1e-12) if dedicated and raw_gb > 0 else 0
    if dedicated:
        sequencing = billed_runs * p["reagent_run_usd"]
        allocated_hours = billed_runs * p["run_hours"]
        billed_gb = billed_runs * p["run_yield_gb"]
    else:
        sequencing = max(p["shared_min_sequencing_usd"], raw_gb * p["shared_sequencing_usd_per_gb"])
        allocated_hours = run_equivalents * p["run_hours"]
        billed_gb = raw_gb
    multiplier = p["assembly_compute_multiplier"] if application == "assembly" else 1
    analysis_hours = config["analysis_setup_hours"] + raw_gb * p["compute_hours_per_gb"] * multiplier
    storage_gb = raw_gb * p["stored_gb_per_sequence_gb"]
    library_rate = p["rna_library_per_sample_usd"] if application == "expression" else p["library_per_sample_usd"]
    components = {
        "library_preparation": n * library_rate,
        "sample_qc": n * config["qc_per_sample_usd"],
        "sequencing_consumables": sequencing,
        "instrument_time": allocated_hours * p["instrument_usd_per_hour"],
        "cloud_compute": analysis_hours * config["compute_hour_usd"],
        "storage": storage_gb * config["storage_gb_month_usd"] * config["storage_months"],
    }
    subtotal = sum(components.values())
    components["contingency"] = subtotal * config["contingency_fraction"]
    total = subtotal + components["contingency"]
    uncertainty = config["uncertainty_fraction"]
    return {
        "currency": "USD", "components_usd": {key: round(value, 2) for key, value in components.items()},
        "total_usd": round(total, 2), "total_unrounded_usd": total,
        "planning_range_usd": [round(total * (1 - uncertainty), 2), round(total * (1 + uncertainty), 2)],
        "range_interpretation": "Sensitivity range, not a statistical confidence interval.",
        "raw_gb": round(raw_gb, 6), "billed_gb": round(billed_gb, 6),
        "flow_cell_equivalents": round(run_equivalents, 6), "dedicated_flow_cells": billed_runs if dedicated else None,
        "allocated_instrument_hours": round(allocated_hours, 4), "analysis_hours": round(analysis_hours, 4),
        "storage_gb": round(storage_gb, 4), "storage_months": config["storage_months"],
        "procurement": scenario["procurement"], "provenance": p["cost_provenance"],
    }


def affordable_raw_gb(scenario: dict, platform: str, config: dict) -> float:
    budget = scenario["budget_usd"]
    if cost_for_raw_gb(0, scenario, platform, config)["total_unrounded_usd"] > budget:
        return 0.0
    lower, upper = 0.0, 1.0
    while cost_for_raw_gb(upper, scenario, platform, config)["total_unrounded_usd"] <= budget:
        upper *= 2
        if upper > 1e12:
            break
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if cost_for_raw_gb(midpoint, scenario, platform, config)["total_unrounded_usd"] <= budget:
            lower = midpoint
        else:
            upper = midpoint
    return lower


def validate_cost_benchmarks(config: dict) -> list[dict]:
    """Compare matching cost scopes; report deviations instead of hiding them.

    Published service amounts include instrument operation. They are compared
    with model sequencing + instrument, never with sequencing reagents alone.
    This is a limited external plausibility check, not observed invoice error.
    """
    validation = []
    for benchmark in config.get("validation_benchmarks", []):
        row = dict(benchmark)
        platform = config["platforms"][row["platform"]]
        kind = row.get("kind", "reagent_arithmetic")
        if kind == "reagent_arithmetic":
            model_unit = platform["reagent_run_usd"] / platform["run_yield_gb"]
            benchmark_unit = row["published_reagent_usd"] / row["raw_gb"]
            row.update(model_reagent_usd_per_gb=round(model_unit, 6),
                       arithmetic_relative_difference_pct=round((model_unit / benchmark_unit - 1) * 100, 4))
        elif kind == "independent_sequencing_service":
            gb_range = [pairs * row["read_bp"] * 2 / 1e9 for pairs in row["read_pairs_range"]]
            unit_range = [row["price_usd"] / gb_range[1], row["price_usd"] / gb_range[0]]
            model_unit = (platform["shared_sequencing_usd_per_gb"] +
                          platform["run_hours"] * platform["instrument_usd_per_hour"] / platform["run_yield_gb"])
            row.update(raw_gb_range=gb_range,
                       service_usd_per_gb_range=[round(value, 6) for value in unit_range],
                       model_sequencing_plus_instrument_usd_per_gb=round(model_unit, 6),
                       relative_difference_pct_range=[round((model_unit / unit_range[1] - 1) * 100, 4),
                                                      round((model_unit / unit_range[0] - 1) * 100, 4)],
                       comparison_status="below_reference_range" if model_unit < unit_range[0] else
                                         "above_reference_range" if model_unit > unit_range[1] else "inside_reference_range")
        elif kind == "independent_preparation_service":
            model_price = platform[row["field"]]
            row.update(model_per_sample_usd=model_price,
                       relative_difference_pct=round((model_price / row["price_usd"] - 1) * 100, 4),
                       comparison_status="below_reference" if model_price < row["price_usd"] else
                                         "above_reference" if model_price > row["price_usd"] else "equal_reference")
        validation.append(row)
    return validation


def make_candidate(scenario: dict, platform: str, config: dict, alt_mean: float) -> dict:
    p = config["platforms"][platform]
    application = scenario["application"]
    retention = scenario["mapping_fraction"] * scenario["usable_fraction"]
    genome_gb = scenario["genome_size_mb"] / 1000
    n = scenario["sample_count"]
    caveats: list[str] = []
    reasons: list[str] = []
    compatible = True
    read_length = p["read_length"]
    read_length_bp = p["read_length_bp"]
    effective_depth: float | None = None
    count_target: float | None = None
    detection: dict[str, Any]
    if application == "variant":
        floor = config["variant_floor_x"][scenario["organism"]]
        target = alt_mean / scenario["allele_fraction"]
        effective_depth = float(math.ceil(max(floor, target)))
        raw_gb = effective_depth * genome_gb * n / retention
        limit = alt_mean / effective_depth
        detection = {
            "type": "idealized_allele_sampling", "target_allele_fraction": scenario["allele_fraction"],
            "minimum_allele_fraction_at_target_probability": min(1, limit),
            "requested_probability": scenario["detection_probability"],
            "predicted_sampling_probability": poisson_tail(effective_depth * scenario["allele_fraction"], scenario["min_alt_reads"]),
            "minimum_alt_reads": scenario["min_alt_reads"], "required_mean_alt_reads": round(alt_mean, 8),
            "specificity_validated": False,
            "interpretation": "Probability of sufficient alternate observations at an accessible locus under ideal independent sampling; not measured caller recall or a trustworthy analytical limit of detection.",
        }
        reasons.append(f"Use the greater of the provisional {floor:g}× design floor and {target:.1f}× required for alternate-read sampling.")
        reasons.append("Short reads are efficient for SNVs in uniquely mappable regions; HiFi can improve access to complex loci.")
        if platform == "nanopore":
            caveats.append("Nanopore systematic errors require chemistry/basecaller-specific truth-set validation, especially for low-frequency variants. The Poisson curve does not model those errors.")
        if scenario["allele_fraction"] < 0.05:
            caveats.append("Rare-allele designs need validated error suppression, unique molecules and contamination controls; ordinary depth alone does not establish specificity.")
    elif application == "assembly":
        effective_depth = float(config["assembly_floor_x"][platform])
        raw_gb = effective_depth * genome_gb * n / retention
        required_span = scenario["repeat_span_kb"] * 1000 + 2000
        if scenario["assembly_goal"] == "resolve_repeats":
            if platform == "illumina" and required_span > p.get("planning_insert_span_bp", 350):
                compatible = False
                reasons.append("This short-read insert cannot bridge the requested repeat plus two 1-kb unique anchors.")
            elif platform == "pacbio" and required_span > read_length_bp:
                compatible = False
                reasons.append("The illustrative 15-kb HiFi library target is shorter than the requested repeat plus anchors.")
            elif platform == "nanopore" and required_span > read_length_bp:
                target_length = p.get("ultralong_target_bp", 100000)
                if required_span > target_length:
                    compatible = False
                    reasons.append("Requested span exceeds the configured ultra-long library target; this model cannot certify a bridging design.")
                else:
                    read_length_bp = target_length
                    read_length = "≥100 kb ultra-long DNA library target (must be confirmed experimentally)"
                    caveats.append("Ultra-long extraction changes yield and cost. The standard DNA cost envelope may underestimate this design.")
        else:
            reasons.append("Short reads can support a draft assembly; repeat-spanning completeness is not requested.")
        detection = {
            "type": "assembly_span_requirement", "required_span_bp": required_span,
            "planning_read_length_bp": read_length_bp, "assembly_completeness_validated": False,
            "interpretation": "Read length plus depth is a planning filter. It does not predict contig N50, phase accuracy or assembly completeness.",
        }
        reasons.append(f"{effective_depth:g}× effective depth is an explicitly provisional assembly target; genome complexity changes it.")
        caveats.append("A single mean or target read length does not guarantee repeat bridging; confirm the full read-length distribution and spanning-read count.")
    else:
        if scenario["expression_goal"] == "gene_counts":
            if platform != "illumina":
                compatible = False
                reasons.append("The implemented gene-count budget model is calibrated in short-read pairs; long-read RNA requires a separate calibrated yield model.")
            count_target = max(config["expression_min_pairs_millions"], scenario["expression_pairs_millions"]) * 1e6
            bases_per_fragment = 300
            read_length = "2 × 150 bp RNA-seq" if platform == "illumina" else "Gene-count short-read model not applicable to this platform"
        else:
            if platform == "illumina":
                compatible = False
                reasons.append("Short reads cannot directly observe complete multi-exon transcript molecules for the requested full-length isoform objective.")
            count_target = max(config["expression_min_long_reads_millions"], scenario["expression_long_reads_millions"]) * 1e6
            bases_per_fragment = 1500
            read_length = "Full-length cDNA; 1.5-kb mean transcript assumption"
            caveats.append("Long-read RNA estimates use an illustrative 1.5-kb mean transcript and DNA-equivalent flow-cell allocation. RNA yield and isoform sensitivity remain unvalidated.")
        raw_gb = count_target * bases_per_fragment / 1e9 * n / retention
        transcript_mean = required_poisson_mean(scenario["detection_probability"], int(config["expression_count_threshold"]))
        detection = {
            "type": "idealized_transcript_sampling", "effective_fragments_per_sample": count_target,
            "minimum_count_threshold": config["expression_count_threshold"],
            "minimum_fragment_fraction_ppm": transcript_mean / count_target * 1e6,
            "requested_probability": scenario["detection_probability"],
            "power_validated": False,
            "interpretation": "Sampling threshold in fraction of assigned fragments, not TPM and not differential-expression statistical power.",
        }
        reasons.append("Expression requires fragment counts and biological replication, not genome coverage; deeper sequencing cannot substitute for independent samples.")
        if scenario["replicates_per_condition"] < 3:
            caveats.append("Fewer than three planned biological replicates per condition is a weak exploratory design; model dispersion and power using pilot counts.")
        if n < scenario["replicates_per_condition"] * scenario["condition_count"]:
            caveats.append("The sample count is smaller than conditions × planned replicates; the biological design is incomplete.")
        caveats.append("Transcript abundance, RNA quality, rRNA depletion, library bias and biological dispersion determine actual sensitivity.")
    cost = cost_for_raw_gb(raw_gb, scenario, platform, config)
    affordable_gb = affordable_raw_gb(scenario, platform, config)
    affordable_depth = affordable_gb * retention / (genome_gb * n) if application != "expression" else None
    affordable_count = None
    if application == "expression":
        affordable_count = affordable_gb * retention * 1e9 / (n * bases_per_fragment)
    over_budget = cost["total_unrounded_usd"] > scenario["budget_usd"] + 1e-8
    adequate = not (application == "expression" and n < scenario["replicates_per_condition"] * scenario["condition_count"])
    return {
        "platform": platform, "label": p["label"], "instrument": p["instrument"],
        "technology": p["technology"], "read_length": read_length, "read_length_bp": read_length_bp,
        "error_mode": p["error_mode"], "capability_note": p["capability_note"], "sources": p["sources"],
        "compatible": compatible, "biological_design_adequate": adequate,
        "feasible": compatible and adequate and not over_budget,
        "within_budget": not over_budget, "budget_gap_usd": round(max(0, cost["total_usd"] - scenario["budget_usd"]), 2),
        "effective_depth_x": effective_depth, "nominal_depth_x": round(effective_depth / retention, 3) if effective_depth is not None else None,
        "effective_fragments_per_sample": count_target,
        "affordable_effective_depth_x": round(affordable_depth, 3) if affordable_depth is not None else None,
        "affordable_fragments_per_sample": round(affordable_count) if affordable_count is not None else None,
        "cost": cost, "detection": detection, "reasons": reasons, "caveats": caveats,
    }


def recommend(scenario: dict[str, Any], config: dict | None = None, evidence: dict | None = None) -> dict:
    normalized = normalize_input(scenario)
    config = load_config() if config is None else config
    evidence = load_evidence() if evidence is None else evidence
    alt_mean = required_poisson_mean(normalized["detection_probability"], normalized["min_alt_reads"])
    platforms = ("illumina", "pacbio", "nanopore") if normalized["platform"] == "auto" else (normalized["platform"],)
    candidates = [make_candidate(normalized, platform, config, alt_mean) for platform in platforms]
    eligible = [candidate for candidate in candidates if candidate["compatible"]]
    feasible = [candidate for candidate in eligible if candidate["feasible"]]
    chosen = min(feasible or eligible or candidates, key=lambda candidate: (
        candidate["cost"]["total_unrounded_usd"], candidate["platform"]))
    if not eligible:
        status = "unsupported"
        headline = "No supported design for the requested objective"
    elif not chosen["biological_design_adequate"]:
        status = "incomplete_design"
        headline = "Complete the biological sample design"
    elif not feasible:
        status = "over_budget"
        headline = "Required design exceeds the planning budget"
    else:
        status = "within_estimate"
        headline = "A design fits the current planning estimate"
    methods = {
        "retention_formula": "effective_fraction = mapping_fraction × usable_fraction",
        "effective_fraction": normalized["mapping_fraction"] * normalized["usable_fraction"],
        "depth_formula": "raw_Gb = effective_depth × haploid_genome_Mb / 1000 × sample_count / effective_fraction",
        "detection_formula": "P(alternate_reads ≥ k) = 1 − Σ(j=0…k−1) exp(−D×f) (D×f)^j / j!",
        "sampling_model": "Poisson depth with independent alternate observations; D is mean mapped, retained depth and f is allele fraction.",
        "cost_formula": "sum(library, QC, sequencing, instrument, compute, storage) × (1 + contingency)",
        "selection": "Lowest estimated total cost among compatible designs meeting the requested budget and sample structure; if none fit, report the least costly compatible requirement and its gap.",
        "floor_basis": "Depth floors and RNA fragment targets are explicit planning priors. The measured bacterial mapping pilot has insufficient endpoint coverage to replace them.",
    }
    detection_curve = []
    if normalized["application"] == "variant":
        target_depth = chosen["effective_depth_x"] or 30
        depths = sorted(set([5, 10, 20, 30, 40, 60, 100, 150, 300, 500, math.ceil(target_depth)]))
        detection_curve = [{"effective_depth_x": depth,
                            "sampling_probability": poisson_tail(depth * normalized["allele_fraction"], normalized["min_alt_reads"])}
                           for depth in depths]
    limitations = [
        "This is a reproducible research planning model. End-to-end variant sensitivity, expression power and assembly performance have not been validated.",
        "Mapping and usable fractions are editable assumptions, not predictions from organism names. Replace them with QC and alignment results from your library.",
        "Mean depth hides uneven coverage, duplicate molecules and inaccessible regions. The allele model excludes sequencing errors, strand bias and caller filters.",
        "Vendor reagent benchmarks check unit arithmetic. Independent NUSeq service checks report configuration/scope differences and substantial preparation-price mismatches; they do not validate the project total.",
        "The cost range is a ±30% sensitivity scenario by default, not a measured error distribution; quotations and atypical sample preparation can fall outside it.",
    ]
    config_bytes = json.dumps(config, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return {
        "schema_version": 1, "model_version": MODEL_VERSION,
        "config_version": config.get("version"), "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "status": status, "headline": headline, "input": normalized, "recommendation": chosen,
        "alternatives": candidates, "detection_curve": detection_curve, "methodology": methods,
        "evidence": evidence, "limitations": limitations,
        "sources": config.get("sources", []), "cost_validation": validate_cost_benchmarks(config),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, help="JSON scenario file")
    parser.add_argument("--output", type=Path, help="Save JSON (default: stdout)")
    parser.add_argument("--cost-config", type=Path, help="Alternative editable cost JSON")
    parser.add_argument("--evidence-dir", type=Path, help="Directory containing results/ JSON files")
    parser.add_argument("--organism", choices=sorted(ENUMS["organism"]))
    parser.add_argument("--application", choices=sorted(ENUMS["application"]))
    parser.add_argument("--budget", type=float)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--genome-mb", type=float)
    parser.add_argument("--allele-fraction", type=float)
    parser.add_argument("--platform", choices=sorted(ENUMS["platform"]))
    parser.add_argument("--procurement", choices=sorted(ENUMS["procurement"]))
    args = parser.parse_args(argv)
    try:
        scenario = read_json(args.input) if args.input else {}
        if not isinstance(scenario, dict):
            raise InputError("Input file must contain a JSON object.")
        overrides = {"organism": args.organism, "application": args.application, "budget_usd": args.budget,
                     "sample_count": args.samples, "genome_size_mb": args.genome_mb,
                     "allele_fraction": args.allele_fraction, "platform": args.platform, "procurement": args.procurement}
        scenario.update({key: value for key, value in overrides.items() if value is not None})
        result = recommend(scenario, load_config(args.cost_config), load_evidence(args.evidence_dir))
        output = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(output, encoding="utf-8")
        else:
            sys.stdout.write(output)
    except (InputError, OSError, ValueError) as error:
        sys.stderr.write(f"Error: {error}\n")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
