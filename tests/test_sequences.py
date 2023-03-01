"""
Template Test Cases, please update
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

class TestTemplate(TestCase):
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
