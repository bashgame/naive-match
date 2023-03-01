"""
Sequence Test Cases
"""
import random
from unittest import TestCase
from sequences import utils

comps = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N'}

def makeSeq(length):
    seq = ''
    for _ in range(length):
        seq += random.choice('ACGT')
    return seq

class TestSequences(TestCase):
    """ Test Template, please update """

    @classmethod
    def setUpClass(cls):
        """ Do something to set up tests, please update """

    @classmethod
    def tearDownClass(cls):
        """ Do something to clean up tests, please update """

    def setUp(self):
        """ Ensure each test runs clean, please update """

    def tearDown(self):
        """ Ensure each test runs clean, please update """

    ###########################################################################
    #   T E S T  C A S E S
    ###########################################################################

    def test_rev_comp(self):
        """ It should return the reverse complement of the sequence """
        test_seq1 = "GATTACA"
        test_seq2 = "gattaca"
        res1 = utils.rev_comp(test_seq1)
        res2 = utils.rev_comp(test_seq2)
        self.assertEqual(res1, 'TGTAATC')
        self.assertEqual(res2, 'TGTAATC')

    def test_naive_match(self):
        """ It should return a list of indices where a complete match occurred """
        test_seq1 = makeSeq(100)
        test_seq2 = test_seq1[1:4]
        expected = test_seq1.rfind(test_seq2)
        result = utils.naive(test_seq2, test_seq1)
        self.assertTrue(expected in result)

    def test_non_match(self):
        """ It should return an empty list """
        test_seq = makeSeq(100)
        bad_seq = "BOB"
        expected = []
        result = utils.naive(bad_seq, test_seq)
        self.assertEqual(expected, result)

    def test_rc_naive(self):
        """ It should return matches for the pattern and its reverse complement """
        test_seq = makeSeq(100)
        test_pat = 'GTA'
        test_rc = utils.rev_comp(test_pat)
        expected = utils.naive(test_pat, test_seq)
        expected += utils.naive(test_rc, test_seq)
        expected = sorted(expected)
        result = sorted(utils.rc_naive(test_pat, test_seq))
        self.assertEqual(expected, result)

    def test_rc_naive_own_rc(self):
        """ It should return an index only once """
        test_seq = makeSeq(1000) + 'AGCT'                   # Ensure at least 1 match
        test_pat = 'AGCT'
        test_rc = utils.rev_comp(test_pat)
        expected = utils.naive(test_pat, test_seq)
        expected += utils.naive(test_rc, test_seq)          # This will double up the matches
        expected = sorted(expected)
        result = sorted(utils.rc_naive(test_pat, test_seq))
        self.assertNotEqual(expected, result)

    def test_get_seq_from_fasta(self):
        """ It should return a sequence """
        filename = 'lambda_virus.fa'
        seq = utils.seq_from_fasta(filename)
        self.assertTrue(len(seq) > 0)
        self.assertFalse(seq.startswith('>'))

    def test_file_not_found(self):
        """ It should return an error message """
        filename = 'dingo_baby'
        errorMessage = 'Error, file not found'
        seq = utils.seq_from_fasta(filename)
        self.assertEqual(seq, errorMessage)

