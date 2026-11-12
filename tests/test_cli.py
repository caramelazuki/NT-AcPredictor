"""CLI regressions. Run: python3 -m unittest discover -s tests -v."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "code/python/nt-AcPredictor.py"
INTRO = (
    "In the Acetylation column, Ac and unAc represents Acetylated and "
    "Unacetylated respectively.\n\n"
)
WARNING = (
    "Since the sequence of which the amino acid in the first site (without iMet) "
    "is not A, C, G, P, S, T, V\n"
    "may not be processed by NatA or B, the prediction result for the sequence "
    "may not be accurate.\n\n"
)
HEADER = (
    "*************** Prediction results ***************\n"
    "ID\tSequence\tAcetylation\tSupplement\n"
)
SHORT = "The following sequence(s) without iMet is shorter than 10 mer.\n"
INVALID = "The following sequence(s) contains non-standard character(s).\n"


class CliTests(unittest.TestCase):
    def run_cli(self, fasta, mode):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fasta"
            path.write_text(fasta, encoding="utf-8")
            return subprocess.run(
                [sys.executable, str(SCRIPT), "-i", str(path), "-f", mode],
                capture_output=True, text=True, timeout=10, check=False,
            )

    def fasta(self, records, mode):
        # "exist" removes the initiator residue before validating the 10-mer.
        prefix = "M" if mode == "exist" else ""
        return "".join(
            ">" + name + "\n" + prefix + sequence + "\n"
            for name, sequence in records
        )

    def assert_input_error(self, records, expected):
        for mode in ("truncated", "exist"):
            with self.subTest(mode=mode):
                result = self.run_cli(self.fasta(records, mode), mode)
                self.assertEqual(result.returncode, 1)
                self.assertEqual(result.stdout, expected)
                self.assertEqual(result.stderr, "")

    def test_malformed_fasta(self):
        for mode in ("truncated", "exist"):
            for fasta in ("AAAAAAAAAAA\n", ">missing_sequence\n"):
                with self.subTest(mode=mode, fasta=fasta):
                    result = self.run_cli(fasta, mode)
                    self.assertEqual(result.returncode, 1)
                    self.assertEqual(
                        result.stdout, "Input file is not of standard FASTA format.\n"
                    )
                    self.assertEqual(result.stderr, "")

    def test_short_sequence(self):
        self.assert_input_error([("short", "AAAAAAAAA")], SHORT + "short\n")

    def test_nonstandard_character(self):
        self.assert_input_error([("invalid", "AAAAAXAAAA")], INVALID + "invalid\n")

    def test_both_errors_preserve_message_order(self):
        self.assert_input_error(
            [("both", "AX"), ("short", "AAA"), ("invalid", "AAAAAXAAAA")],
            SHORT + "both\nshort\n" + INVALID + "both\ninvalid\n",
        )

    def test_invalid_record_prevents_partial_predictions(self):
        self.assert_input_error(
            [("valid", "AAAAAAAAAA"), ("invalid", "AAAAAXAAAA")],
            INVALID + "invalid\n",
        )

    def test_valid_predictions(self):
        # Fixed outputs verified against the pre-fix program, including exact
        # 10-residue inputs, longer/wrapped sequences, and duplicate 10-mers.
        records = [
            ("acetylated", "AAAAAAAAAA"),
            ("unacetylated", "APAAAAAAAA"),
            ("wrapped", "AAAAA\nAAAAACDE"),
            ("duplicate", "AAAAAAAAAA"),
        ]
        expected = INTRO + HEADER + (
            "acetylated\tAAAAAAAAAA\tAc\t\n"
            "unacetylated\tAPAAAAAAAA\tunAc\t\n"
            "wrapped\tAAAAAAAAAA\tAc\t\n"
            "duplicate\tAAAAAAAAAA\tAc\t\n"
        )
        for mode in ("truncated", "exist"):
            with self.subTest(mode=mode):
                result = self.run_cli(self.fasta(records, mode), mode)
                self.assertEqual(result.returncode, 0)
                self.assertEqual(result.stdout, expected)
                self.assertEqual(result.stderr, "")

    def test_valid_prediction_warning(self):
        expected = INTRO + WARNING + HEADER + (
            "warning\tRAAAAAAAAA\tunAc\t"
            "The prediction may not be accurate. The sequence may be processed "
            "by not NatA or D but NatB, C E or F.\n"
        )
        for mode in ("truncated", "exist"):
            with self.subTest(mode=mode):
                result = self.run_cli(
                    self.fasta([("warning", "RAAAAAAAAA")], mode), mode
                )
                self.assertEqual(result.returncode, 0)
                self.assertEqual(result.stdout, expected)
                self.assertEqual(result.stderr, "")

    def test_help(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True, timeout=10, check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("-i INPUT", result.stdout)
        self.assertIn("-f {truncated,exist}", result.stdout)
        self.assertIn("yamada@tohoku.ac.jp", result.stdout)
        self.assertEqual(result.stderr, "")


if __name__ == "__main__":
    unittest.main()
