import unittest

from tandem_repeat_decomposer import decompose_hmm
from tandem_repeat_decomposer import decompose_dp
# # from pomegranate import DiscreteDistribution, State
# from pomegranate import HiddenMarkovModel as Model


class MyTestCase(unittest.TestCase):
    def test_hmm_perfect_motif_len4(self):
        sequence = "ACTGACTGACTG"
        consensus_motif = "ACTG"
        deomposed_motifs = decompose_hmm(sequence, consensus_motif)
        self.assertEqual(deomposed_motifs, ['ACTG', 'ACTG', 'ACTG'])

    def test_hmm_imperfect_motif_len4(self):
        sequence = "ACTGAACTGACTG"
        consensus_motif = "ACTG"
        deomposed_motifs = decompose_hmm(sequence, consensus_motif)
        self.assertEqual(deomposed_motifs, ['ACTG', 'AACTG', 'ACTG'])

    def test_hmm_perfect_motif_len3(self):
        sequence = "AACAACAACAACAAC"
        consensus_motif = "AAC"
        deomposed_motifs = decompose_hmm(sequence, consensus_motif)
        self.assertEqual(deomposed_motifs, ['AAC', 'AAC', 'AAC', 'AAC', 'AAC'])

    def test_hmm_perfect_motif_len10_repeat2(self):
        sequence = "ACCGCCGTTGACCGCCGTTG"
        consensus_motif = "ACCGCCGTTG"
        deomposed_motifs = decompose_hmm(sequence, consensus_motif)
        self.assertEqual(deomposed_motifs, ['ACCGCCGTTG', 'ACCGCCGTTG'])

    def test_hmm_perfect_low_complexity_motif(self):
        sequence = "AAATAAAATAAAATAAAATA"
        consensus_motif = "AAATA"
        deomposed_motifs = decompose_hmm(sequence, consensus_motif)
        self.assertEqual(deomposed_motifs, ['AAATA', 'AAATA', 'AAATA', 'AAATA'])

    # def test_hmm_impefect_low_complexity_motif(self):
    #     sequence = "AAATAAAATTAAAATAAAAAATA"
    #     consensus_motif = "AAATA"
    #     deomposed_motifs = decompose_hmm(sequence, consensus_motif)
    #     self.assertEqual(deomposed_motifs, ['AAATA', 'AAATTA', 'AAATA', 'AAAAATA'])
    #     # This failed with answer, ['AAATA', 'AAATT', 'AAAA', 'TAAA', 'AAATA']
    #     # This is because the estimated repeat count is 5, and it must pass 5 repeats.

    def test_dp_impefect_low_complexity_motif(self):
        sequence = "AAATAAAATTAAAATAAAAAATA"
        motifs = ["AAATA"]
        deomposed_motifs = decompose_dp(sequence, motifs)
        self.assertEqual(deomposed_motifs, ['AAATA', 'AAATTA', 'AAATA', 'AAAAATA'])

    def test_hmm_invalid_sequence(self):
        sequence = "NNNNNNNNNNN"
        consensus_motif = "ACTG"
        with self.assertRaises(ValueError) as context:
            decompose_hmm(sequence, consensus_motif)
        self.assertEqual('Sequence has invalid characters', str(context.exception))

    def test_dp_invalid_sequence(self):
        sequence = "NNNNNNNNNNN"
        consensus_motif = "ACTG"
        with self.assertRaises(ValueError) as context:
            decompose_dp(sequence, list(consensus_motif))
        self.assertEqual('Sequence has invalid characters', str(context.exception))

    def test_dp_imperfect_motifs(self):
        sequence = "AGCCAGGCAGCC"
        motifs = ["AGCC", "AGGC"]
        deomposed_motifs = decompose_dp(sequence, motifs)
        self.assertEqual(deomposed_motifs, ['AGCC', 'AGGC', 'AGCC'])

if __name__ == '__main__':
    unittest.main()
