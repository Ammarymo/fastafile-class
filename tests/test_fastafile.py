import unittest
from fastafile.fastafile import FastaFile

class TestFastaFile(unittest.TestCase):
    def setUp(self):
        self.fasta = FastaFile('examples/example.fasta')

    def test_get_entries(self):
        entries = self.fasta.get_entries()
        self.assertEqual(len(entries), 2)

    def test_calculate_gc_content(self):
        sequence = "ATCGATCG"
        gc_content = self.fasta.calculate_gc_content(sequence)
        self.assertAlmostEqual(gc_content, 50.0)

    def test_reverse_complement(self):
        sequence = "ATCG"
        reverse_complement = self.fasta.reverse_complement(sequence)
        self.assertEqual(reverse_complement, "CGAT")

if __name__ == "__main__":
    unittest.main()
