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
