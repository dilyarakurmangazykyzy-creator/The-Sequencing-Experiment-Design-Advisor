"""Synthetic fixtures test correctness only; never substitute for real evidence."""
import gzip
import importlib.util
from pathlib import Path
import random
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pipeline as p


class FastqTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def write(self, name, text):
        path = self.directory / name
        with gzip.open(path, "wt") as handle:
            handle.write(text)
        return path

    def test_qc_and_gzip_stream(self):
        path = self.write("reads.fastq.gz", "@r1\nACGTNN\n+\nIIIIII\n@r2\nAAAA\n+\n!!!!\n")
        qc = p.QC()
        for read in p.fastq_records(path):
            qc.add(read)
        report = qc.report()
        self.assertEqual(report["reads"], 2)
        self.assertEqual(report["bases"], 10)
        self.assertAlmostEqual(report["gc_fraction_acgt"], .25)
        self.assertAlmostEqual(report["n_base_fraction"], .2)
        self.assertAlmostEqual(report["q30_base_fraction"], .6)
        self.assertEqual(report["read_length_median"], 4)

    def test_truncated_fastq_rejected(self):
        path = self.write("bad.fastq.gz", "@r1\nACGT\n+\n")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            list(p.fastq_records(path))

    def test_mismatched_quality_rejected(self):
        path = self.write("bad.fastq.gz", "@r1\nACGT\n+\nIII\n")
        with self.assertRaisesRegex(ValueError, "length mismatch"):
            list(p.fastq_records(path))

    def test_pair_ids_and_count(self):
        first = self.write("r1.fastq.gz", "@a/1\nACGT\n+\nIIII\n")
        second = self.write("r2.fastq.gz", "@a/2\nACGT\n+\nIIII\n")
        self.assertEqual(len(list(p.read_units(first, second))), 1)
        wrong = self.write("wrong.fastq.gz", "@b/2\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "disagree"):
            list(p.read_units(first, wrong))
        empty = self.write("empty.fastq.gz", "")
        with self.assertRaisesRegex(ValueError, "different numbers"):
            list(p.read_units(first, empty))

    def test_platform_specific_filtering(self):
        read = p.Read("r", "A" * 60 + p.ADAPTER + "A" * 30, "I" * 103)
        filtered, reason = p.filter_read(read, "illumina")
        self.assertIsNone(reason)
        self.assertEqual(len(filtered.sequence), 60)
        long = p.Read("n", "A" * 600, "(" * 600)  # Q7
        self.assertIsNotNone(p.filter_read(long, "nanopore")[0])
        low = p.Read("n", "A" * 600, "'" * 600)  # Q6
        self.assertEqual(p.filter_read(low, "nanopore")[1], "mean_phred_lt_7")


class SamplingTests(unittest.TestCase):
    def test_pair_selection_and_seed_reproducibility(self):
        for i in range(200):
            self.assertEqual(p.selected(f"read{i}/1", .4, 17), p.selected(f"read{i}/2", .4, 17))
            self.assertEqual(p.hash_fraction(f"read{i}", 17), p.hash_fraction(f"read{i}", 17))
        ids = [f"read{i}" for i in range(1000)]
        a = {i for i in ids if p.selected(i, .1, 17)}
        b = {i for i in ids if p.selected(i, .5, 17)}
        self.assertTrue(a <= b)
        self.assertNotEqual(a, {i for i in ids if p.selected(i, .1, 29)})
        self.assertTrue(all(p.selected(i, 1, 17) for i in ids))
        self.assertFalse(any(p.selected(i, 0, 17) for i in ids))


class MetadataTests(unittest.TestCase):
    def test_missing_metadata_and_read_denominator(self):
        spec = importlib.util.spec_from_file_location("metadata_mine", ROOT / "scripts" / "metadata_mine.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            path.write_text("run_accession\tinstrument_platform\tlibrary_strategy\tlibrary_layout\tfirst_public\tread_count\tbase_count\n"
                            "SRR1\tILLUMINA\tWGS\tPAIRED\t2025-01-01\t20\t2000\n"
                            "SRR2\tILLUMINA\tWGS\tPAIRED\t2025-01-01\t\t90000\n", encoding="utf-8")
            result = module.aggregate(path, directory, quiet=True)
            self.assertEqual(result["rows_processed"], 2)
            self.assertEqual(result["missing_counts"]["reads"], 1)
            self.assertEqual(result["missing_counts"]["spots"], 2)
            self.assertEqual(result["design_groups"][0]["aggregate_bases_per_archive_count_unit"], 100)
            self.assertEqual(result["design_groups"][0]["spots"], 0)
            self.assertIn("spots", result["schema_unavailable_fields"])
            self.assertNotIn("spots", result["blank_counts_given_returned_column"])
            self.assertEqual(result["blank_counts_given_returned_column"]["reads"], 1)


class AlignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(122)
        cls.reference = "".join(rng.choice("ACGT") for _ in range(5000))
        cls.mapper = p.SeedMapper([cls.reference], step=1)

    def test_unique_reverse_and_unmapped(self):
        fragment = self.reference[1200:1350]
        hit = self.mapper.map_short(fragment)
        self.assertTrue(hit["mapped"])
        self.assertTrue(hit["unique"])
        self.assertEqual(hit["intervals"], [(1200, 1350)])
        self.assertIsNone(hit["mapq"])
        reverse = self.mapper.map_short(p.reverse_complement(fragment))
        self.assertEqual(reverse["strand"], "-")
        self.assertEqual(reverse["intervals"], [(1200, 1350)])
        self.assertFalse(self.mapper.map_short("N" * 150)["mapped"])

    def test_repeat_ambiguity_excluded_from_breadth(self):
        fragment = self.reference[1100:1350]
        repeated = p.SeedMapper([self.reference[:300] + fragment + self.reference[300:600] + fragment], step=1)
        result = repeated.map_short(fragment[30:180])
        self.assertTrue(result["mapped"])
        self.assertFalse(result["unique"])
        self.assertEqual(result["intervals"], [])

    def test_short_substitution_verification(self):
        fragment = self.reference[500:650]
        replacement = "A" if fragment[30] != "A" else "C"
        edited = fragment[:30] + replacement + fragment[31:]
        result = self.mapper.map_short(edited)
        self.assertTrue(result["mapped"])
        self.assertAlmostEqual(result["error_fraction"], 1 / 150)

    def test_banded_edit_correctness(self):
        self.assertEqual(p.banded_edit_distance("ACGT", "ACGT", 1), 0)
        self.assertEqual(p.banded_edit_distance("ACGT", "AGGT", 1), 1)
        self.assertEqual(p.banded_edit_distance("ACGT", "ACGGT", 1), 1)
        self.assertGreater(p.banded_edit_distance("AAAA", "TTTT", 2), 2)

    def test_long_exact_anchor_breadth_is_conservative(self):
        fragment = self.reference[1800:3000]
        result = self.mapper.map_long(fragment)
        self.assertTrue(result["mapped"])
        covered = sum(e - s for s, e in result["intervals"])
        self.assertGreater(covered, 500)
        self.assertLess(covered, len(fragment))
        self.assertIsNone(result["error_fraction"])

    def test_sam_primary_cigar_and_mapq(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mapping.sam"
            path.write_text("@SQ\tSN:chr\tLN:100\n"
                            "r1\t0\tchr\t11\t40\t5M2I3M2D4M\t*\t0\t0\tACGTACGTACGTAA\tIIIIIIIIIIIIII\tNM:i:4\n"
                            "r2\t0\tchr\t21\t5\t4M\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n"
                            "r3\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n"
                            "unknownq\t0\tchr\t11\t255\t4M\t*\t0\t0\tAAAA\tIIII\n"
                            "duplicate\t1024\tchr\t11\t40\t4M\t*\t0\t0\tAAAA\tIIII\n"
                            "r1\t256\tchr\t40\t30\t4M\t*\t0\t0\tAAAA\tIIII\n", encoding="ascii")
            rows = p.read_sam(path, {"chr": 0}, {"chr": 100})
        self.assertEqual(len(rows["r1"]), 1)
        covered = {pos for start, end in rows["r1"][0]["intervals"] for pos in range(start, end)}
        self.assertEqual(covered, set(range(10, 18)) | set(range(20, 24)))
        self.assertEqual(rows["r2"][0]["intervals"], [])
        self.assertTrue(rows["r2"][0]["mapped"])
        self.assertFalse(rows["r3"][0]["mapped"])
        self.assertEqual(rows["unknownq"][0]["intervals"], [])
        self.assertIsNone(rows["unknownq"][0]["mapq"])
        self.assertEqual(rows["duplicate"][0]["intervals"], [])
        self.assertAlmostEqual(rows["r1"][0]["error_fraction"], 4 / 14)

    def test_reference_breadth_and_depth(self):
        records = [{"intervals": [(0, 4), (2, 6)]}, {"intervals": [(2, 3)]}]
        report = p.coverage_metrics(records, 10)
        self.assertAlmostEqual(report["coverage_1x"], .6)
        self.assertEqual(report["coverage_5x"], 0)
        self.assertAlmostEqual(report["mapped_mean_depth_x"], .9)
        original_np = p.np
        try:
            p.np = None
            fallback = p.coverage_metrics(records, 10)
        finally:
            p.np = original_np
        self.assertAlmostEqual(fallback["coverage_1x"], report["coverage_1x"])
        self.assertAlmostEqual(fallback["mapped_mean_depth_x"], report["mapped_mean_depth_x"])

    def test_gc_windows_boundary_dropout_and_invalid_coordinates(self):
        records = [{"intervals": [(0, 4), (2, 6)]}, {"intervals": [(2, 3)]}]
        original_np = p.np
        try:
            p.np = None
            result = p.coverage_metrics(records, 10, [(0, 5, .2), (5, 10, .8)])
        finally:
            p.np = original_np
        self.assertAlmostEqual(result["gc_dropout"][0]["mean_depth"], 1.6)
        self.assertAlmostEqual(result["gc_dropout"][0]["uncovered_fraction"], 0)
        self.assertAlmostEqual(result["gc_dropout"][1]["mean_depth"], .2)
        self.assertAlmostEqual(result["gc_dropout"][1]["uncovered_fraction"], .8)
        if original_np is not None:
            vector = p.coverage_metrics(records, 10, [(0, 5, .2), (5, 10, .8)])
            self.assertEqual(vector["gc_dropout"], result["gc_dropout"])
        with self.assertRaisesRegex(ValueError, "outside reference"):
            p.coverage_metrics([{"intervals": [(-1, 4)]}], 10)
        with self.assertRaisesRegex(ValueError, "GC window"):
            p.coverage_metrics(records, 10, [(5, 11, .5)])
        with self.assertRaisesRegex(ValueError, "disjoint"):
            p.coverage_metrics(records, 10, [(0, 6, .5), (5, 10, .5)])

    def test_sam_contig_boundary_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mapping.sam"
            path.write_text("r\t0\tchr\t98\t40\t4M\t*\t0\t0\tAAAA\tIIII\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "beyond reference"):
                p.read_sam(path, {"chr": 0}, {"chr": 100})

    def test_quality_means_handle_unavailable_unmapped_and_empty_reads(self):
        rows = [{"mapped": True, "mapq": 40, "error_fraction": .01},
                {"mapped": True, "mapq": 0, "error_fraction": .05},
                {"mapped": True, "mapq": 255, "error_fraction": None},
                {"mapped": False, "mapq": 60, "error_fraction": .5}]
        scores = p.alignment_quality_metrics(rows)
        self.assertEqual(scores["mean_mapq"], 20)
        self.assertAlmostEqual(scores["mean_alignment_edit_distance_per_aligned_query_base"], .03)
        self.assertEqual(scores["scored_mapq_reads"], 2)
        self.assertEqual(scores["scored_edit_distance_reads"], 2)
        empty = p.alignment_quality_metrics([])
        self.assertIsNone(empty["mean_mapq"])
        self.assertIsNone(empty["mean_alignment_edit_distance_per_aligned_query_base"])
        self.assertEqual(empty["scored_mapq_reads"], 0)


if __name__ == "__main__":
    unittest.main()
