"""Model contract and scientific-unit tests: python -m unittest discover -s tests."""

from __future__ import annotations

import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import advisor


class AdvisorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = advisor.load_config()
        cls.no_evidence = {"status": "absent", "saturation": None, "summary": None}

    def recommend(self, **input):
        return advisor.recommend(input, self.config, self.no_evidence)

    def test_invalid_numeric_inputs(self):
        invalid = [
            {"budget_usd": 0}, {"budget_usd": -1}, {"sample_count": 0},
            {"sample_count": 1.5}, {"genome_size_mb": 0}, {"allele_fraction": 0},
            {"allele_fraction": 1.1}, {"mapping_fraction": 0}, {"usable_fraction": 0},
            {"detection_probability": 1}, {"min_alt_reads": 0}, {"budget_usd": True},
            {"sample_count": "6"}, {"unexpected": 3}, {"platform": "magic"},
        ]
        for scenario in invalid:
            with self.subTest(scenario=scenario):
                with self.assertRaises(advisor.InputError):
                    self.recommend(**scenario)

    def test_nan_and_infinity_rejected_in_every_numeric_field(self):
        for field in advisor.NUMBERS:
            for value in (float("nan"), float("inf"), -float("inf")):
                with self.subTest(field=field, value=value):
                    with self.assertRaises(advisor.InputError):
                        self.recommend(**{field: value})

    def test_poisson_detection_reference_value(self):
        mean = advisor.required_poisson_mean(0.95, 3)
        self.assertAlmostEqual(mean, 6.29579362187, places=9)
        self.assertAlmostEqual(advisor.poisson_tail(mean, 3), 0.95, places=12)
        self.assertAlmostEqual(advisor.required_poisson_mean(0.95, 1), -math.log(0.05), places=10)
        self.assertEqual(advisor.poisson_tail(0, 3), 0)
        self.assertAlmostEqual(advisor.poisson_tail(10000, 100), 1)

    def test_low_allele_fraction_needs_more_depth(self):
        high = self.recommend(organism="human", sample_count=1, allele_fraction=0.5)["recommendation"]
        low = self.recommend(organism="human", sample_count=1, allele_fraction=0.05)["recommendation"]
        self.assertGreater(low["effective_depth_x"], high["effective_depth_x"])
        self.assertEqual(low["effective_depth_x"], 126)
        self.assertGreaterEqual(low["detection"]["predicted_sampling_probability"], 0.95)
        self.assertFalse(low["detection"]["specificity_validated"])

    def test_budget_infeasible_has_explicit_gap(self):
        result = self.recommend(budget_usd=1)
        self.assertEqual(result["status"], "over_budget")
        self.assertFalse(result["recommendation"]["feasible"])
        self.assertGreater(result["recommendation"]["budget_gap_usd"], 0)
        self.assertEqual(result["recommendation"]["affordable_effective_depth_x"], 0)

    def test_budget_monotonicity_and_requirement_independence(self):
        low = self.recommend(platform="illumina", budget_usd=800)["recommendation"]
        high = self.recommend(platform="illumina", budget_usd=1600)["recommendation"]
        self.assertGreaterEqual(high["affordable_effective_depth_x"], low["affordable_effective_depth_x"])
        self.assertEqual(high["effective_depth_x"], low["effective_depth_x"])
        self.assertEqual(high["cost"]["total_usd"], low["cost"]["total_usd"])

    def test_dedicated_runs_round_up_with_unused_capacity(self):
        input = advisor.normalize_input({"sample_count": 1, "procurement": "dedicated"})
        capacity = self.config["platforms"]["nanopore"]["run_yield_gb"]
        for raw, expected in ((capacity * 0.5, 1), (capacity, 1), (capacity + 0.01, 2)):
            cost = advisor.cost_for_raw_gb(raw, input, "nanopore", self.config)
            self.assertEqual(cost["dedicated_flow_cells"], expected)
            self.assertEqual(cost["billed_gb"], expected * capacity)
            self.assertEqual(cost["components_usd"]["sequencing_consumables"], expected * 840)

    def test_nominal_and_effective_depth_units(self):
        result = self.recommend(sample_count=2, genome_size_mb=5, mapping_fraction=0.8, usable_fraction=0.5)
        design = result["recommendation"]
        self.assertAlmostEqual(design["nominal_depth_x"], design["effective_depth_x"] / 0.4)
        expected_gb = design["effective_depth_x"] * 0.005 * 2 / 0.4
        self.assertAlmostEqual(design["cost"]["raw_gb"], expected_gb, places=6)

    def test_paired_read_length_counts_both_mates(self):
        result = self.recommend(application="expression", sample_count=6, expression_pairs_millions=30,
                                mapping_fraction=1, usable_fraction=1)
        self.assertEqual(result["recommendation"]["platform"], "illumina")
        # 30 million pairs × (150+150) bp × 6 = 54 Gb, not 27 Gb.
        self.assertAlmostEqual(result["recommendation"]["cost"]["raw_gb"], 54)
        self.assertIsNone(result["recommendation"]["effective_depth_x"])
        different_genome = self.recommend(application="expression", genome_size_mb=4.6,
                                         sample_count=6, expression_pairs_millions=30,
                                         mapping_fraction=1, usable_fraction=1)
        self.assertEqual(result["recommendation"]["cost"]["raw_gb"],
                         different_genome["recommendation"]["cost"]["raw_gb"])

    def test_sample_replication_and_repeat_span_constraints(self):
        expression = self.recommend(application="expression", sample_count=2)
        self.assertEqual(expression["status"], "incomplete_design")
        assembly = self.recommend(application="assembly", platform="illumina", repeat_span_kb=8)
        self.assertEqual(assembly["status"], "unsupported")
        draft = self.recommend(application="assembly", platform="illumina", assembly_goal="draft", budget_usd=10000)
        self.assertTrue(draft["recommendation"]["compatible"])

    def test_independent_service_scope_and_pair_units(self):
        rows = advisor.validate_cost_benchmarks(self.config)
        lane = next(row for row in rows if row.get("name", "").startswith("NUSeq 10B"))
        self.assertEqual(lane["raw_gb_range"], [300, 360])
        self.assertAlmostEqual(lane["service_usd_per_gb_range"][0], 1873.75 / 360, places=6)
        # Include model instrument allocation when comparing with a service price.
        self.assertEqual(lane["model_sequencing_plus_instrument_usd_per_gb"], 5.15)
        self.assertEqual(lane["comparison_status"], "below_reference_range")
        prep = next(row for row in rows if row.get("name", "").startswith("NUSeq DNA"))
        self.assertAlmostEqual(prep["relative_difference_pct"], -66.6667, places=4)

    def test_missing_and_unrecognized_evidence_not_fabricated(self):
        with tempfile.TemporaryDirectory() as directory:
            empty = advisor.load_evidence(directory)
            self.assertEqual(empty["status"], "absent")
            output = Path(directory) / "results"
            output.mkdir()
            (output / "saturation.json").write_text(json.dumps({"endpoint":"snv_recall", "platforms":{}}))
            rejected = advisor.load_evidence(directory)
            self.assertEqual(rejected["status"], "absent")
            self.assertTrue(any("Unrecognized" in text for text in rejected["issues"]))

    def test_deterministic_result(self):
        first = self.recommend(allele_fraction=0.1)
        second = self.recommend(allele_fraction=0.1)
        self.assertEqual(first, second)
        self.assertEqual(len(first["config_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
