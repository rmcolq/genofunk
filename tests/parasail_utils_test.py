import os
import unittest
import json
from Bio.Seq import Seq
import filecmp

from genofunk.parasail_utils import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class TestParasailUtils(unittest.TestCase):
    def test_pairwise_sw_align(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "cccatgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataac"
        result = pairwise_sw_align(ref_seq, query_seq)
        (ref_begin, ref_end, read_begin, read_end) = get_alignment_start_end(result)
        self.assertEqual(3,read_begin)
        self.assertEqual(read_end, 56)
        print(result)

    def test_pairwise_nw_trace_align(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        result = pairwise_nw_trace_align(ref_seq, query_seq)
        self.assertEqual(b'8=1X7=1X6=1D9=1X10=1I10=',result.cigar.decode)

    def test_parse_cigar(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        result = pairwise_nw_trace_align(ref_seq, query_seq)
        parsed_cigar = parse_cigar_pairs(result)
        expected = [("=",8), ("X",1), ("=",7), ("X",1), ("=",6), ("D",1), ("=",9), ("X",1), ("=",10), ("I",1), ("=",10)]
        self.assertEqual(expected, parsed_cigar)

    def test_get_cigar_length(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        result = pairwise_nw_trace_align(ref_seq, query_seq)
        pairs = parse_cigar_pairs(result)
        cigar_length = get_cigar_length(pairs, max_mismatch=3)
        expected = 23
        self.assertEqual(expected, cigar_length)

    def test_cigar_length_max_mismatch(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaatacggtagagggttttcaacatttgaggaccgatgtataa"
        result = pairwise_nw_trace_align(ref_seq, query_seq)
        pairs = parse_cigar_pairs(result)
        print(pairs)
        cigar_length = get_cigar_length(pairs, max_mismatch=1)
        expected = 16
        self.assertEqual(expected, cigar_length)

    def test_cigar_length_n_runs(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaaccNNNNNaccgtagagggttttcaacatttgaggaccgatgtataa"
        result = pairwise_nw_trace_align(ref_seq, query_seq)
        pairs = parse_cigar_pairs(result)
        cigar_length = get_cigar_length(pairs, max_mismatch=3, n_runs=[[10,14]], min_match=1)
        expected = 23
        self.assertEqual(expected, cigar_length)

    def test_is_extended_cigar_prefix(self):
        c1 = [("=", 12), ("X", 1), ("=", 17)]
        c2 = [("=", 10), ("X", 3), ("=", 17)]
        c3 = [("=", 12), ("X", 3), ("=", 15)]
        c4 = [("=", 12), ("I", 1), ("=", 17)]
        c5 = [("=", 12), ("D", 1), ("=", 17)]
        c6 = [("=", 12)]

        self.assertEqual(is_extended_cigar_prefix(c1, c1), False)
        self.assertEqual(is_extended_cigar_prefix(c2, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c3, c3), False)
        self.assertEqual(is_extended_cigar_prefix(c4, c4), False)
        self.assertEqual(is_extended_cigar_prefix(c5, c5), False)
        self.assertEqual(is_extended_cigar_prefix(c6, c6), False)

        self.assertEqual(is_extended_cigar_prefix(c6, c1), True)
        self.assertEqual(is_extended_cigar_prefix(c6, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c6, c3), True)
        self.assertEqual(is_extended_cigar_prefix(c6, c4), True)
        self.assertEqual(is_extended_cigar_prefix(c6, c5), True)

        self.assertEqual(is_extended_cigar_prefix(c1, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c2, c1), False)
        self.assertEqual(is_extended_cigar_prefix(c1, c3), False)
        self.assertEqual(is_extended_cigar_prefix(c3, c1), False)
        self.assertEqual(is_extended_cigar_prefix(c1, c4), False)
        self.assertEqual(is_extended_cigar_prefix(c4, c1), False)
        self.assertEqual(is_extended_cigar_prefix(c1, c5), False)
        self.assertEqual(is_extended_cigar_prefix(c5, c1), False)
        self.assertEqual(is_extended_cigar_prefix(c2, c3), False)
        self.assertEqual(is_extended_cigar_prefix(c3, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c2, c4), False)
        self.assertEqual(is_extended_cigar_prefix(c4, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c2, c5), False)
        self.assertEqual(is_extended_cigar_prefix(c5, c2), False)
        self.assertEqual(is_extended_cigar_prefix(c3, c4), False)
        self.assertEqual(is_extended_cigar_prefix(c4, c3), False)
        self.assertEqual(is_extended_cigar_prefix(c3, c5), False)
        self.assertEqual(is_extended_cigar_prefix(c5, c3), False)
        self.assertEqual(is_extended_cigar_prefix(c4, c5), False)
        self.assertEqual(is_extended_cigar_prefix(c5, c4), False)

    def test_is_improved_cigar_prefix(self):
        c1 = [("=", 12), ("X", 1), ("=", 17)]
        c2 = [("=", 10), ("X", 3), ("=", 17)]
        c3 = [("=", 12), ("X", 3), ("=", 15)]
        c4 = [("=", 12), ("I", 1), ("=", 17)]
        c5 = [("=", 12), ("D", 1), ("=", 17)]
        c6 = [("=", 12)]

        self.assertEqual(is_improved_cigar_prefix(c1, c1), False)
        self.assertEqual(is_improved_cigar_prefix(c2, c2), False)
        self.assertEqual(is_improved_cigar_prefix(c3, c3), False)
        self.assertEqual(is_improved_cigar_prefix(c4, c4), False)
        self.assertEqual(is_improved_cigar_prefix(c5, c5), False)
        self.assertEqual(is_improved_cigar_prefix(c6, c6), False)

        self.assertEqual(is_improved_cigar_prefix(c1, c2), False)
        self.assertEqual(is_improved_cigar_prefix(c2, c1), True)
        self.assertEqual(is_improved_cigar_prefix(c1, c3), False)
        self.assertEqual(is_improved_cigar_prefix(c3, c1), True)
        self.assertEqual(is_improved_cigar_prefix(c1, c4), False)
        self.assertEqual(is_improved_cigar_prefix(c4, c1), True)
        self.assertEqual(is_improved_cigar_prefix(c1, c5), False)
        self.assertEqual(is_improved_cigar_prefix(c5, c1), True)
        self.assertEqual(is_improved_cigar_prefix(c6, c1), True)
        self.assertEqual(is_improved_cigar_prefix(c1, c6), False)

        self.assertEqual(is_improved_cigar_prefix(c3, c2), False)
        self.assertEqual(is_improved_cigar_prefix(c2, c3), True)
        self.assertEqual(is_improved_cigar_prefix(c4, c2, True), False)
        self.assertEqual(is_improved_cigar_prefix(c4, c2, False), False)
        self.assertEqual(is_improved_cigar_prefix(c2, c4), True)
        self.assertEqual(is_improved_cigar_prefix(c5, c2), False)
        self.assertEqual(is_improved_cigar_prefix(c2, c5), True)
        self.assertEqual(is_improved_cigar_prefix(c6, c2), False) # funny case, prefix is better, overall is worse
        self.assertEqual(is_improved_cigar_prefix(c2, c6), True)

        self.assertEqual(is_improved_cigar_prefix(c3, c4), False)
        self.assertEqual(is_improved_cigar_prefix(c4, c3), True)
        self.assertEqual(is_improved_cigar_prefix(c3, c5), False)
        self.assertEqual(is_improved_cigar_prefix(c5, c3), True)
        self.assertEqual(is_improved_cigar_prefix(c3, c6), False)
        self.assertEqual(is_improved_cigar_prefix(c6, c3), True)

        self.assertEqual(is_improved_cigar_prefix(c4, c5), False)
        self.assertEqual(is_improved_cigar_prefix(c5, c4), False)
        self.assertEqual(is_improved_cigar_prefix(c4, c6), False)
        self.assertEqual(is_improved_cigar_prefix(c6, c4), True)

        self.assertEqual(is_improved_cigar_prefix(c6, c5), True)
        self.assertEqual(is_improved_cigar_prefix(c5, c6), False)

    def test_is_improved_cigar_prefix_real_examples(self):
        c1 = [('=', 5), ('X', 1), ('=', 5), ('X', 19)]
        c2 = [('=', 5), ('X', 1), ('=', 4), ('X', 1), ('=', 9), ('X', 1), ('=', 9)]
        self.assertEqual(is_improved_cigar_prefix(c1, c2), True)

        c3 = [('=', 102), ('X', 1), ('=', 122)]
        c4 = [('=', 102), ('X', 1), ('=', 121), ('X', 1)]
        self.assertEqual(is_improved_cigar_prefix(c3, c4), False)
