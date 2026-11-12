import contextlib
import io
import unittest

from genome_compare.motif import (
    DEFAULT_MOTIF,
    MotifHit,
    describe,
    find_motif,
    main,
    reverse_complement,
)

INSTANCE = "TATACCCATGGGGT"  # matches TATANNNATGNNNT


class ReverseComplementTests(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(reverse_complement("AACG"), "CGTT")

    def test_lowercase(self):
        self.assertEqual(reverse_complement("acgt"), "ACGT")

    def test_is_its_own_inverse(self):
        self.assertEqual(reverse_complement(reverse_complement(INSTANCE)), INSTANCE)

    def test_n_is_kept(self):
        self.assertEqual(reverse_complement("ACNT"), "ANGT")

    def test_invalid_characters_raise(self):
        with self.assertRaises(ValueError):
            reverse_complement("ACXT")


class FindMotifTests(unittest.TestCase):
    def test_forward_strand_hit(self):
        hits = find_motif("GGCC" + INSTANCE + "AACC")
        self.assertEqual(hits, [MotifHit("+", 5, 18, INSTANCE)])

    def test_reverse_strand_hit_uses_forward_coordinates(self):
        sequence = "GGCC" + reverse_complement(INSTANCE) + "AACC"
        hits = find_motif(sequence)
        self.assertEqual(hits, [MotifHit("-", 5, 18, INSTANCE)])

    def test_no_hit_on_either_strand(self):
        self.assertEqual(find_motif("GGCCTATACCCATGCCCAAAAGGTT"), [])

    def test_lowercase_input(self):
        self.assertEqual(len(find_motif(("GGCC" + INSTANCE).lower())), 1)

    def test_n_in_sequence_does_not_match_a_wildcard(self):
        # NNN in the sequence is not "any three DNA bases"
        self.assertEqual(find_motif("TATANNNATGCCCT"), [])

    def test_overlapping_hits_are_found(self):
        # "AAA" occurs at positions 1 and 2 of "AAAA", so overlaps must count.
        hits = find_motif("AAAA", motif="AAA")
        forward = [h for h in hits if h.strand == "+"]
        self.assertEqual([h.start for h in forward], [1, 2])

    def test_custom_motif_with_wildcard(self):
        hits = find_motif("CCGAATTCCC", motif="GANNTC")
        self.assertIn(MotifHit("+", 3, 8, "GAATTC"), hits)

    def test_invalid_motif_raises(self):
        with self.assertRaises(ValueError):
            find_motif("ACGT", motif="AXG")
        with self.assertRaises(ValueError):
            find_motif("ACGT", motif="")

    def test_invalid_sequence_raises(self):
        with self.assertRaises(ValueError):
            find_motif("ACGTX")

    def test_hits_are_sorted_by_position(self):
        sequence = reverse_complement(INSTANCE) + "CC" + INSTANCE
        starts = [h.start for h in find_motif(sequence)]
        self.assertEqual(starts, sorted(starts))
        self.assertEqual(len(starts), 2)


class DescribeAndCliTests(unittest.TestCase):
    def test_describe_reports_not_found(self):
        self.assertIn("not found", describe("GGCCTATACCCATGCCCAAAAGGTT"))

    def test_describe_reports_strand_and_positions(self):
        text = describe("GGCC" + INSTANCE + "AACC")
        self.assertIn("forward strand", text)
        self.assertIn("5-18", text)

    def test_cli_with_sequence(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = main(["GGCC" + INSTANCE + "AACC"])
        self.assertEqual(code, 0)
        self.assertIn(DEFAULT_MOTIF, out.getvalue())

    def test_cli_invalid_sequence_returns_error(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = main(["ACGTX"])
        self.assertEqual(code, 1)
        self.assertIn("error:", err.getvalue())


if __name__ == "__main__":
    unittest.main()
