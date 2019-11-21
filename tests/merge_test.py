import os
import json
import unittest
from unittest import mock

import filecmp

from genofunk import merge
from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'merge')

class TestMerge(unittest.TestCase):
    def setUp(self):
        self.m = merge.Merge()

    def test_load_coordinates_file_does_not_exist(self):
        coordinate_file = os.path.join(data_dir, "doesnt_exist.coordinates")
        self.assertRaises(FileNotFoundError, self.m.load_coordinates_file, coordinate_file)

    def test_load_coordinates_file_empty(self):
        coordinate_file = os.path.join(data_dir, "empty.coordinates")
        self.assertRaises(json.JSONDecodeError, self.m.load_coordinates_file, coordinate_file)

    def test_load_coordinates_file_no_feature_list(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = None
        self.m.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "ORF2": {"seq2": {"start": 0, "end": 126}, "seq3": {"start": 10, "end": 136}}}
        self.assertEqual(self.m.coordinates, expect_coordinates)

    def test_load_coordinates_file_with_feature_list(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = ["ORF1", "idontexistinasample"]
        self.m.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        self.assertEqual(self.m.coordinates, expect_coordinates)

    def test_load_from_directory_empty(self):
        empty_dir = os.path.join(data_dir, "empty")
        self.assertRaises(AssertionError, self.m.load_from_directory, empty_dir)

    def test_load_from_directory_missing_pair(self):
        missing_dir = os.path.join(data_dir, "missing_consensus")
        self.assertRaises(AssertionError, self.m.load_from_directory, missing_dir)

    def test_load_from_directory(self):
        self.m.load_from_directory(data_dir)
        self.assertIsNotNone(self.m.edits)
        self.assertIsNotNone(self.m.consensus_sequence)
        self.assertEqual(len(self.m.edits.edits), 14)
        self.assertEqual(len(self.m.consensus_sequence), 12)

    def test_query_edit(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        e = self.m.edits.edits[0]
        print(e)
        with mock.patch('builtins.input', return_value="y"):
            self.m.query_edit(e)
            assert e.edit_accepted == False

        with mock.patch('builtins.input', return_value="n"):
            self.m.query_edit(e)
            assert e.edit_accepted == True