import filecmp
import os
import unittest
import glob
import json
from Bio.Seq import Seq

from genofunk import annotator

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'annotator')

class TestAnnotator(unittest.TestCase):
    def setUp(self):
        ref_filepath = os.path.join(data_dir, 'ref.json')
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        self.a = annotator.Annotator("LN854563.1")
        self.a.load_input_files(ref_filepath, consensus_filepath)

    def test_load_reference_info_no_file(self):
        ref_filepath = os.path.join(data_dir, 'idontexist.json')
        a = annotator.Annotator("accession")
        self.assertRaises(FileNotFoundError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_empty_file(self):
        ref_filepath = os.path.join(data_dir, 'empty_ref.json')
        a = annotator.Annotator("accession")
        self.assertRaises(json.decoder.JSONDecodeError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_references(self):
        ref_filepath = os.path.join(data_dir, 'no_references_ref.json')
        a = annotator.Annotator("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_sequence(self):
        ref_filepath = os.path.join(data_dir, 'missing_sequence_ref.json')
        a = annotator.Annotator("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_orf(self):
        ref_filepath = os.path.join(data_dir, 'missing_orf_ref.json')
        a = annotator.Annotator("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_missing_accession(self):
        ref_filepath = os.path.join(data_dir, 'missing_accession_ref.json')
        a = annotator.Annotator("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_correct_accession(self):
        ref_filepath = os.path.join(data_dir, 'ref.json')
        a = annotator.Annotator("LN854563.1")
        data = a.load_reference_info(ref_filepath)
        print(a)
        print(a.closest_accession)
        print(data)
        self.assertIsNotNone(data)

    def test_load_consensus_no_file(self):
        consensus_filepath = os.path.join(data_dir, 'idontexist.fasta')
        a = annotator.Annotator("accession")
        self.assertRaises(FileNotFoundError, a.load_consensus_sequence, consensus_filepath)

    def test_load_consensus_empty_file(self):
        consensus_filepath = os.path.join(data_dir, 'empty_consensus.fasta')
        a = annotator.Annotator("accession")
        self.assertRaises(AssertionError, a.load_consensus_sequence, consensus_filepath)

    def test_load_consensus_simple_case(self):
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        a = annotator.Annotator("accession")
        records = a.load_consensus_sequence(consensus_filepath)
        self.assertEqual(len(records), 8)
        self.assertEqual(records[0].seq, "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaa"
                                          "gatgcgcatattacctaa")
        self.assertEqual(records[1].seq,"aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaag"
                                         "attggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")

    def test_load_input_files(self):

        self.assertIsNotNone(self.a.reference_info)
        self.assertIsNotNone(self.a.consensus_sequence)
        self.assertIsNotNone(self.a.edits)

    def test_get_sequence_not_Seq(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotator.Annotator("LN854563.1")
        self.assertRaises(TypeError, a.get_sequence, seq)

    def test_get_sequence_simple_case(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, amino_acid=False)
        self.assertEqual(seq, result)

    def test_get_sequence_coordinates_na(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, coordinates=(0,12), amino_acid=False)
        expected = "atgcccaagctg"
        self.assertEqual(expected, result)

    def test_get_sequence_offset_na(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, offset=3, amino_acid=False)
        expected = "cccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        self.assertEqual(expected, result)

    def test_get_sequence_aa(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEqual(expected, result)

    def test_get_sequence_coordinates_aa(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, coordinates=(0,48))
        expected = "MPKLNSVEGFSSFEDD"
        self.assertEqual(expected, result)

    def test_get_sequence_aa_length_remainder_1(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtat")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDVX"
        self.assertEqual(expected, result)

    def test_get_sequence_aa_length_remainder_2(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtata")
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDVX"
        self.assertEqual(expected, result)

    def test_get_reference_sequence_no_accession(self):
        result = self.a.get_reference_sequence()
        expected ="MINAHLEINTHEGRNDTHERELIVEDAHIT*NTANASTYDIRTYWETHLEFILLEDWITHTHEENDSFWRMS*"
        self.assertEqual(expected, result)

    def test_get_reference_sequence_accession(self):
        result = self.a.get_reference_sequence(accession="test")
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEqual(expected, result)

    def test_get_query_sequence_no_id(self):
        result = self.a.get_query_sequence(amino_acid=False)
        expected = "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa"
        self.assertEqual(expected, result)

    def test_get_query_sequence_id_1(self):
        result = self.a.get_query_sequence(record_id=1, amino_acid=False)
        expected = "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccat" \
                   "gaagaaaacgatagcttttggcgcatgagctaa"
        self.assertEqual(expected, result)

    # No tests for decode_cigar

    def test_pairwise_ssw_align(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "cccatgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataac"
        a = annotator.Annotator()
        result = a.pairwise_ssw_align(ref_seq, query_seq)
        self.assertEqual(3,result.read_begin1)
        self.assertEqual(result.read_end1, 56)
        print(result)

    def test_pairwise_sw_trace_align(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        a = annotator.Annotator()
        result = a.pairwise_sw_trace_align(ref_seq, query_seq)
        self.assertEqual(b'8=1X7=1X6=1D9=1X10=1I10=',result.cigar.decode)

    def test_parse_cigar(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        a = annotator.Annotator()
        result = a.pairwise_sw_trace_align(ref_seq, query_seq)
        parsed_cigar = a.parse_cigar(result)
        expected = [("=",8), ("X",1), ("=",7), ("X",1), ("=",6), ("D",1), ("=",9), ("X",1), ("=",10), ("I",1), ("=",10)]
        self.assertEqual(expected, parsed_cigar)

    def test_cigar_length(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "atgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataa"
        a = annotator.Annotator()
        result = a.pairwise_sw_trace_align(ref_seq, query_seq)
        pairs = a.parse_cigar(result)
        cigar_length = a.cigar_length(pairs)
        expected = 23
        self.assertEqual(expected, cigar_length)

    def test_is_improved_cigar_prefix(self):
        a = annotator.Annotator()
        c1 = [("=", 12), ("X", 1), ("=", 17)]
        c2 = [("=", 10), ("X", 3), ("=", 17)]
        c3 = [("=", 12), ("X", 3), ("=", 15)]
        c4 = [("=", 12), ("I", 1), ("=", 17)]
        c5 = [("=", 12), ("D", 1), ("=", 17)]
        c6 = [("=", 12)]

        self.assertEqual(a.is_improved_cigar_prefix(c1, c1), False)
        self.assertEqual(a.is_improved_cigar_prefix(c2, c2), False)
        self.assertEqual(a.is_improved_cigar_prefix(c3, c3), False)
        self.assertEqual(a.is_improved_cigar_prefix(c4, c4), False)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c5), False)
        self.assertEqual(a.is_improved_cigar_prefix(c6, c6), False)

        self.assertEqual(a.is_improved_cigar_prefix(c1,c2), False)
        self.assertEqual(a.is_improved_cigar_prefix(c2, c1), True)

        self.assertEqual(a.is_improved_cigar_prefix(c1, c3), False)
        self.assertEqual(a.is_improved_cigar_prefix(c3, c1), True)

        self.assertEqual(a.is_improved_cigar_prefix(c1, c4), False)
        self.assertEqual(a.is_improved_cigar_prefix(c4, c1), True)

        self.assertEqual(a.is_improved_cigar_prefix(c1, c5), False)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c1), True)

        self.assertEqual(a.is_improved_cigar_prefix(c1, c5), False)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c1), True)

        self.assertEqual(a.is_improved_cigar_prefix(c6, c1), True)
        self.assertEqual(a.is_improved_cigar_prefix(c1, c6), False)

        self.assertEqual(a.is_improved_cigar_prefix(c4, c2), False)
        self.assertEqual(a.is_improved_cigar_prefix(c2, c4), True)

        self.assertEqual(a.is_improved_cigar_prefix(c5, c2), False)
        self.assertEqual(a.is_improved_cigar_prefix(c2, c5), True)

        self.assertEqual(a.is_improved_cigar_prefix(c6, c2), False) # funny case, prefix is better, overall is worse
        self.assertEqual(a.is_improved_cigar_prefix(c2, c6), True)

        self.assertEqual(a.is_improved_cigar_prefix(c3, c4), False)
        self.assertEqual(a.is_improved_cigar_prefix(c4, c3), True)

        self.assertEqual(a.is_improved_cigar_prefix(c3, c5), False)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c3), True)

        self.assertEqual(a.is_improved_cigar_prefix(c6, c5), True)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c6), False)

        self.assertEqual(a.is_improved_cigar_prefix(c4, c5), False)
        self.assertEqual(a.is_improved_cigar_prefix(c5, c4), False)

    def test_identify_orf_coordinates(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "cccatgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataac"
        a = annotator.Annotator()
        result = a.pairwise_ssw_align(ref_seq, query_seq)
        self.assertEqual(3, result.read_begin1)
        self.assertEqual(result.read_end1, 56)

    def test_frame_shift_insert_n_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 2
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 12)]
        shift_from = ""
        shift_to = "N"
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1]+1)
        self.assertEqual(result_cigar_pairs, [("=", 12), ("X", 1), ("=", 17)])
        self.assertEqual(result_updated, True)

    def test_frame_shift_insert_n_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 2
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 14)]
        shift_from = ""
        shift_to = "N"
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1])
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_insert_nn_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 3
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = self.a.pairwise_sw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = self.a.parse_cigar(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 20)]
        shift_from = ""
        shift_to = "NN"
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1]+2)
        self.assertEqual(result_cigar_pairs, [("=", 20), ("X", 1), ("=", 9)])
        self.assertEqual(result_updated, True)

    def test_frame_shift_insert_nn_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 3
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 21)]
        shift_from = ""
        shift_to = "NN"
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1])
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_delete_n_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 4
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = self.a.pairwise_sw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = self.a.parse_cigar(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 7)]
        shift_from = "N"
        shift_to = ""
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1]-1)
        self.assertEqual(self.a.cigar_length(result_cigar_pairs),29)
        self.assertEqual(result_updated, True)

    def test_frame_shift_delete_n_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 4
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 9)]
        shift_from = "N"
        shift_to = ""
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1])
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_delete_nn_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 5
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = self.a.pairwise_sw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = self.a.parse_cigar(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 12)]
        shift_from = "NN"
        shift_to = ""
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1]-2)
        self.assertEqual(self.a.cigar_length(result_cigar_pairs),29)
        self.assertEqual(result_updated, True)

    def test_frame_shift_delete_nn_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = 5
        ref_sequence = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 13)]
        shift_from = "NN"
        shift_to = ""
        (result_found_coordinates, result_cigar_pairs, result_updated) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record_id,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        shift_from,
                                                                                        shift_to)
        self.assertEqual(result_found_coordinates[0], found_coordinates[0])
        self.assertEqual(result_found_coordinates[1], found_coordinates[1])
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_discover_edits_mismatches_no_edits(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        record_id = 6
        self.a.discover_edits(orf_coordinates, found_coordinates, record_id)

        print(self.a.edits.edits)
        self.assertEqual(len(self.a.edits.edits),0)

    def test_discover_edits_mismatch_insertion_deletion_mismatch(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        record_id = 7
        self.a.discover_edits(orf_coordinates, found_coordinates, record_id)

        print(self.a.edits.edits)
        self.assertEqual(len(self.a.edits.edits),2)