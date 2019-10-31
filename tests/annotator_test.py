import filecmp
import os
import unittest
import glob
import json

from genofunk import annotator

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'annotator')

class TestAnnotator(unittest.TestCase):
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
        ref_filepath = os.path.join(data_dir, 'missing_accession_ref.json')
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
        self.assertEquals(len(records), 2)
        self.assertEquals(records[0].seq, "ATAAGCCCATAAAAATTGAGAACACCCACATGGCCCTTAACTTAATTAGTAGGGGTCCATCTCCAATTCCACAAGATAAACCCCCTAAAGACCAGAGGGATAAACCCCCAAGGAATGT")
        self.assertEquals(records[1].seq,"AGGCAATGGGATGGATCGACCCCCCCATGGACCAGAATCTACCAACCTGGGAGGAGTTAAGCCAAACAGAAAAGCAAGAGATACTCAAAAACAATT")

    def test_load_input_files(self):
        ref_filepath = os.path.join(data_dir, 'missing_accession_ref.json')
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        a = annotator.Annotator("LN854563.1")
        a.load_input_files(ref_filepath, consensus_filepath)
        self.assertIsNotNone(a.reference_info)
        self.assertIsNotNone(a.consensus_sequence)
        self.assertIsNotNone(a.edits)

    def get_sequence_simple_case(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, amino_acid=False)
        self.assertEquals(seq, result)

    def get_sequence_coordinates(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq, )
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEquals(expected, result)

    def get_sequence_translated(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEquals(expected, result)

    def get_sequence_coordinates(self):
        seq = "atgcccaagctgaatagcgtagaggggttttcatcatttgaggacgatgtataa"
        a = annotator.Annotator("LN854563.1")
        result = a.get_sequence(seq)
        expected = "MPKLNSVEGFSSFEDDV*"
        self.assertEquals(expected, result)