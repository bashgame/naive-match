"""
Template Test Cases, please update
"""
from unittest import TestCase
from sequences import utils

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
