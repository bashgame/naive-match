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

    def test_naive_2mm(self):
        """ It should allow 2 mismatches per occurrence """
        test_seq = makeSeq(100)
        test_pat = test_seq[:4]
        expected = test_seq.find(test_pat)
        test_pat2 = test_pat[0] + utils.rev_comp(test_pat[1]) + utils.rev_comp(test_pat[2]) + test_pat[3]
        test_pat3 = test_pat[0] + utils.rev_comp(test_pat[1]) + utils.rev_comp(test_pat[2]) + utils.rev_comp(test_pat[3])
        result2 = utils.naive_2mm(test_pat2, test_seq)
        result3 = utils.naive_2mm(test_pat3, test_seq)
        self.assertTrue(expected in result2)
        self.assertTrue(expected not in result3)

    def test_get_reads_from_fastq(self):
        """ It should return a dictionary of reads """
        filename = 'ERR037900_1.first1000.fastq'
        reads, quals = utils.reads_from_fastq(filename)
        self.assertTrue(len(reads) > 0 and len(quals) > 0)

    def test_file_not_found(self):
        """ It should return an error message """
        filename = 'dingo_baby'
        errorMessage = 'Error, file not found'
        seq = utils.reads_from_fastq(filename)
        self.assertEqual(seq, errorMessage)

    def test_phred33_to_Q(self):
        """ It should return the numerical equivalent of the phred33 character """
        test_char = '#'
        expected = 2
        result = utils.phred33toQ(test_char)
        self.assertEqual(expected, result)