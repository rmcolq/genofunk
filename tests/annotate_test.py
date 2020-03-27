import os
import unittest
import json
from Bio.Seq import Seq
import filecmp

from genofunk import annotate
from genofunk import editfile
from genofunk.parasail_utils import *
from genofunk.sequence_utils import *

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'annotate')


class TestAnnotate(unittest.TestCase):
    def setUp(self):
        ref_filepath = os.path.join(data_dir, 'ref.json')
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        self.a = annotate.Annotate("hobbit")
        self.a.load_input_files(ref_filepath, consensus_filepath)

    def test_load_reference_info_no_file(self):
        ref_filepath = os.path.join(data_dir, 'idontexist.json')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_empty_file(self):
        ref_filepath = os.path.join(data_dir, 'empty_ref.json')
        a = annotate.Annotate("accession")
        self.assertRaises(json.decoder.JSONDecodeError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_references(self):
        ref_filepath = os.path.join(data_dir, 'no_references_ref.json')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_features(self):
        ref_filepath = os.path.join(data_dir, 'no_features_ref.json')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_features_empty(self):
        ref_filepath = os.path.join(data_dir, 'features_empty_ref.json')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_missing_accession(self):
        ref_filepath = os.path.join(data_dir, 'ref.json')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_sequence(self):
        ref_filepath = os.path.join(data_dir, 'missing_sequence_ref.json')
        a = annotate.Annotate("hobbit")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_no_locations(self):
        ref_filepath = os.path.join(data_dir, 'missing_locations_ref.json')
        a = annotate.Annotate("hobbit")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_locations_empty(self):
        ref_filepath = os.path.join(data_dir, 'locations_empty_ref.json')
        a = annotate.Annotate("test")
        self.assertRaises(AssertionError, a.load_reference_info, ref_filepath)

    def test_load_reference_info_correct_accession(self):
        ref_filepath = os.path.join(data_dir, 'ref.json')
        a = annotate.Annotate("hobbit")
        data = a.load_reference_info(ref_filepath)
        self.assertIsNotNone(data)

    def test_load_consensus_no_file(self):
        consensus_filepath = os.path.join(data_dir, 'idontexist.fasta')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_consensus_sequence, consensus_filepath)

    def test_load_consensus_empty_file(self):
        consensus_filepath = os.path.join(data_dir, 'empty_consensus.fasta')
        a = annotate.Annotate("accession")
        self.assertRaises(AssertionError, a.load_consensus_sequence, consensus_filepath)

    def test_load_consensus_simple_case(self):
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        a = annotate.Annotate("accession")
        a.load_consensus_sequence(consensus_filepath)
        records = [a.consensus_sequence[record_id] for record_id in a.consensus_sequence]
        self.assertEqual(len(records), 10)
        self.assertEqual(records[0].seq, "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaa"
                                          "gatgcgcatattacctaa")
        self.assertEqual(records[1].seq,"aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaag"
                                         "attggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        a.consensus_sequence.close()

    def test_load_input_files(self):
        self.assertIsNotNone(self.a.reference_info)
        self.assertIsNotNone(self.a.consensus_sequence)
        self.assertIsNotNone(self.a.edits)

    def test_get_sequence_not_Seq(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotate.Annotate("hobbit")
        self.assertRaises(TypeError, get_sequence, seq)

    def test_get_sequence_simple_case(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq, amino_acid=False)
        self.assertEqual(seq, result)

    def test_get_sequence_coordinates_na(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq, coordinates=(0,12), amino_acid=False)
        expected = "atgcccaagctg"
        self.assertEqual(expected, result)

    def test_get_sequence_offset_na(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq, offset=3, amino_acid=False)
        expected = "cccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        self.assertEqual(expected, result)

    def test_get_sequence_aa(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEqual(expected, result)

    def test_get_sequence_coordinates_aa(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq, coordinates=(0,48))
        expected = "MPKLNSVEGFSSFEDD"
        self.assertEqual(expected, result)

    def test_get_sequence_aa_length_remainder_1(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtat")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDVX"
        self.assertEqual(expected, result)

    def test_get_sequence_aa_length_remainder_2(self):
        seq = Seq("atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtata")
        a = annotate.Annotate("hobbit")
        result, coordinates = get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDVX"
        self.assertEqual(expected, result)

    def test_get_reference_sequence_no_accession(self):
        result, coordinates = self.a.get_reference_sequence()
        expected ="MINAHLEINTHEGRNDTHERELIVEDAHIT*NTANASTYDIRTYWETHLEFILLEDWITHTHEENDSFWRMS*"
        self.assertEqual(expected, result)

    def test_get_reference_sequence_accession(self):
        result, coordinates = self.a.get_reference_sequence(accession="test")
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEqual(expected, result)

    def test_get_reference_sequence_join(self):
        json_value = self.a.reference_info["references"]["test2"]["locations"]["ORF1"]
        self.a.closest_accession = "test2"
        coordinates = get_coordinates_from_json(json_value, pairs=False)
        result, coordinates = self.a.get_reference_sequence(accession="test2", coordinates=coordinates, amino_acid=False)
        expected = "atgcccaagcttgaatagcgt"
        self.assertEqual(expected, result)

    def test_get_query_sequence_id_1(self):
        result, coordinates = self.a.get_query_sequence("seq2", amino_acid=False)
        expected = "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccat" \
                   "gaagaaaacgatagcttttggcgcatgagctaa"
        self.assertEqual(expected, result)

    # No tests for decode_cigar

    def test_identify_orf_coordinates(self):
        ref_seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        query_seq = "cccatgcccaacctgaataccgtagagggttttcaacatttgaggaccgatgtataac"
        a = annotate.Annotate()
        result = pairwise_sw_align(ref_seq, query_seq)
        (ref_begin, ref_end, read_begin, read_end) = get_alignment_start_end(result)
        self.assertEqual(3, read_begin)
        self.assertEqual(read_end, 56)

    def test_get_position_for_frame_shift_no_indel_or_stop(self):
        found_coordinates = (0, 54)
        record_id = "seq1"
        stop_codons = ["*"]
        max_mismatch = 3
        include_compensatory = False

        query_seq = self.a.consensus_sequence[record_id].seq
        query_aa, c = get_sequence(query_seq, amino_acid=True)

        ref_seq =   Seq("attaacgcgcatcGggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(position)
        self.assertEqual(position, 30)

        ref_seq = Seq("attaacgcgcatGtggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(position)
        self.assertEqual(position, 30)

        ref_seq = Seq("attaacgcgcatctCgaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(position)
        self.assertEqual(position, 30)

    def test_get_position_for_frame_shift_insertion(self):
        found_coordinates = (0, 54)
        record_id = "seq1"
        stop_codons = ["*"]
        max_mismatch = 3
        include_compensatory = False

        query_seq = self.a.consensus_sequence[record_id].seq
        query_aa, c = get_sequence(query_seq, amino_acid=True)

#                      attaacgcgCatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 3)

#                      attaacgcgcAtctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgctctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 3)
#                      attaacgcgcatctggaaattaacacccaTgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcatctggaaattaacacccagaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 9)

    def test_get_position_for_frame_shift_insertion2(self):
        found_coordinates = (0, 54)
        record_id = "seq1"
        stop_codons = ["*"]
        max_mismatch = 3
        include_compensatory = False

        query_seq = self.a.consensus_sequence[record_id].seq
        query_aa, c = get_sequence(query_seq, amino_acid=True)

#                      attaacgcgcatCtggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcattggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 4)

#                      attaacgcgcatcTggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcatcggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 4)
#                      attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 30)

    def test_get_position_for_frame_shift_deletion(self):
        found_coordinates = (0, 54)
        record_id = "seq1"
        stop_codons = ["*"]
        max_mismatch = 3
        include_compensatory = False

        query_seq = self.a.consensus_sequence[record_id].seq
        query_aa, c = get_sequence(query_seq, amino_acid=True)

#                        attaacgcgcat ctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq =   Seq("attaacgcgcattctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 4)


#                      attaacgcgcatc tggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcatcttggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 5) #synonymous

#                      attaacgcgcatct ggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa
        ref_seq = Seq("attaacgcgcatcttggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        ref_aa, c = get_sequence(ref_seq, amino_acid=True)
        result = pairwise_nw_trace_align(ref_aa, query_aa)
        cigar_pairs = parse_cigar_pairs(result)
        position = self.a.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons,
                                                       max_mismatch)
        print(query_aa)
        print(ref_aa)
        print(cigar_pairs)
        print(position)
        self.assertEqual(position, 5) #synonymous

    def test_frame_shift_insert_n_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_1_deletion"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 12)]
        stop_codons = ['*']
        shift_from = ""
        shift_to = "N"
        shift_position = 12
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 1)
        self.assertEqual(result_cigar_pairs, [("=", 12), ("X", 1), ("=", 17)])
        self.assertEqual(result_updated, True)

    def test_frame_shift_insert_n_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_1_deletion"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 14)]
        stop_codons = ['*']
        shift_from = ""
        shift_to = "N"
        shift_position = 14
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 0)
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_insert_nn_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_2_deletions"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence, coordinates = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = parse_cigar_pairs(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 20)]
        stop_codons = ['*']
        shift_from = ""
        shift_to = "NN"
        shift_position = 20
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 2)
        self.assertEqual(result_cigar_pairs, [("=", 20), ("X", 1), ("=", 9)])
        self.assertEqual(result_updated, True)

    def test_frame_shift_insert_nn_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_2_deletions"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 22)]
        stop_codons = ['*']
        shift_from = ""
        shift_to = "NN"
        shift_position = 22
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 0)
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_delete_n_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_1_insertion"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence, coordinates = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = parse_cigar_pairs(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 7)]
        stop_codons = ['*']
        shift_from = "N"
        shift_to = ""
        shift_position = 7
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, -1)
        self.assertEqual(get_position_first_indel_or_mismatch_in_cigar(result_cigar_pairs, max_mismatch=3), 29)
        self.assertEqual(result_updated, True)

    def test_frame_shift_delete_n_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_1_insertion"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 9)]
        stop_codons = ['*']
        shift_from = "N"
        shift_to = ""
        shift_position = 9
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 0)
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_frame_shift_delete_nn_sticks(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_2_insertions"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)

        query_sequence, coordinates = self.a.get_query_sequence(record_id, coordinates=found_coordinates)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = parse_cigar_pairs(result)
        print(cigar_pairs)

        cigar_pairs = [("=", 12)]
        stop_codons = ['*']
        shift_from = "NN"
        shift_to = ""
        shift_position = 12
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, -2)
        self.assertEqual(get_position_first_indel_or_mismatch_in_cigar(result_cigar_pairs, max_mismatch=3), 29)
        self.assertEqual(result_updated, True)

    def test_frame_shift_delete_nn_rejected_worse(self):
        orf_coordinates = (3,93)
        found_coordinates = (0,90)
        record_id = "seq1_with_2_insertions"
        record = self.a.consensus_sequence[record_id]
        ref_sequence, coordinates = self.a.get_reference_sequence(coordinates=orf_coordinates)
        cigar_pairs = [("=", 14)]
        stop_codons = ['*']
        shift_from = "NN"
        shift_to = ""
        shift_position = 14
        (result_coordinate_difference, result_cigar_pairs, result_updated, edit) = self.a.frame_shift(orf_coordinates,
                                                                                        found_coordinates,
                                                                                        record,
                                                                                        ref_sequence,
                                                                                        cigar_pairs,
                                                                                        stop_codons,
                                                                                        shift_from,
                                                                                        shift_to,
                                                                                        shift_position,
                                                                                        include_compensatory=False)
        self.assertEqual(result_coordinate_difference, 0)
        self.assertEqual(result_cigar_pairs, cigar_pairs)
        self.assertEqual(result_updated, False)

    def test_discover_frame_shift_edits_mismatches_no_edits(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        stop_codons = ["*"]
        max_mismatch = 3
        include_compensatory = False
        record_id = "seq1_with_mismatches"
        record = self.a.consensus_sequence[record_id]
        self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, max_mismatch, include_compensatory, record_id)
        self.assertEqual(len(self.a.edits.edits),0)

    def test_discover_frame_shift_edits_mismatch_insertion_deletion_mismatch_include_compensatory(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        stop_codons = ["*"]
        max_mismatch = 1
        include_compensatory = True
        record_id = "seq1_with_mismatch_insertion_deletion_mismatch"
        record = self.a.consensus_sequence[record_id]
        self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, max_mismatch, include_compensatory, record_id)
        self.assertEqual(len(self.a.edits.edits),2)

    def test_discover_frame_shift_edits_mismatch_insertion_deletion_mismatch(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        stop_codons = ["*"]
        max_mismatch = 1
        include_compensatory = False
        record_id = "seq1_with_mismatch_insertion_deletion_mismatch"
        record = self.a.consensus_sequence[record_id]
        self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, max_mismatch, include_compensatory, record_id)
        self.assertEqual(len(self.a.edits.edits),2) # doesn't have to find these, but does

    def test_discover_frame_shift_edits_double_deletion_mismatch_insertion_include_compensatory(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        stop_codons = ["*"]
        max_mismatch = 1
        include_compensatory = True
        record_id = "seq1_with_double_deletion_mismatch_insertion"
        record = self.a.consensus_sequence[record_id]
        self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, max_mismatch, include_compensatory, record_id)
        self.a.edits.sort()
        self.assertEqual(len(self.a.edits.edits),2)

    def test_discover_frame_shift_edits_double_deletion_mismatch_insertion(self):
        orf_coordinates = (3, 93)
        found_coordinates = (0, 90)
        stop_codons = ["*"]
        max_mismatch = 1
        include_compensatory = False
        record_id = "seq1_with_double_deletion_mismatch_insertion"
        record = self.a.consensus_sequence[record_id]
        self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, max_mismatch, include_compensatory, record_id)
        self.a.edits.sort()
        self.assertEqual(len(self.a.edits.edits),2) # doesn't have to find these, but does

    # def test_discover_frame_shift_edits_3_insertions(self):
    #     orf_coordinates = (93, 219)
    #     found_coordinates = (0, 126)
    #     stop_codons = ["*"]
    #     record_id = 9
    #     record_id = "seq2_with_3_insertions"
    #     record = self.a.consensus_sequence[record_id]
    #     self.a.discover_frame_shift_edits(orf_coordinates, found_coordinates, stop_codons, record_id)
    #     self.a.edits.sort()
    #     self.assertEqual(len(self.a.edits.edits),3)
    #     self.assertEqual(self.a.edits.edits[0].reference_position, 115)
    #
    # def test_run(self):
    #     a = annotate.Annotate("hobbit")
    #     ref_filepath = os.path.join(data_dir, 'ref.json')
    #     consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
    #     print("running")
    #     a.run(ref_filepath, consensus_filepath)
    #     print("finished running")
    #
    #     tmp_file = os.path.join(data_dir, 'consensus.fasta.edits')
    #     expect_file = os.path.join(data_dir, 'expect_consensus.fasta.edits')
    #     self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
    #     os.unlink(tmp_file)
    #
    #     tmp_file = os.path.join(data_dir, 'consensus.fasta.coordinates')
    #     expect_file = os.path.join(data_dir, 'expect_consensus.fasta.coordinates')
    #     self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
    #     os.unlink(tmp_file)

    def tearDown(self):
        self.a.consensus_sequence.close()
