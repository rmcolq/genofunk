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
        self.m.coordinates = {}
        self.m.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "ORF2": {"seq2": {"start": 0, "end": 126}, "seq3": {"start": 10, "end": 136}}}
        self.assertEqual(self.m.coordinates, expect_coordinates)

    def test_load_coordinates_file_with_feature_list(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = ["ORF1", "idontexistinasample"]
        self.m.coordinates = {}
        self.m.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        self.assertEqual(self.m.coordinates, expect_coordinates)

    def test_load_edits_in_range_no_coordinates(self):
        edit_file = os.path.join(data_dir, "simple_edits")
        self.assertRaises(AssertionError, self.m.load_edits_in_range, edit_file)

    def test_load_edits_in_range_no_feature_list(self):
        edit_file = os.path.join(data_dir, "simple_edits")
        features_list = None
        self.m.edits = EditFile()
        self.m.coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        self.m.load_edits_in_range(edit_file, features_list)
        expect_edits = EditFile(edit_file)
        self.assertEqual(self.m.edits, expect_edits)

    def test_load_edits_in_range_with_feature_list(self):
        edit_file = os.path.join(data_dir, "simple_edits")
        features_list = ["ORF1", "idontexistinasample"]
        self.m.edits = EditFile()
        self.m.coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        self.m.load_edits_in_range(edit_file, features_list)
        expect_edit_file = os.path.join(data_dir, "simple_edits_expect")
        expect_edits = EditFile(expect_edit_file)
        self.assertEqual(self.m.edits, expect_edits)

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
        self.assertIsNotNone(self.m.coordinates)
        self.assertEqual(len(self.m.edits.edits), 14)
        self.assertEqual(len(self.m.consensus_sequence), 12)
        print(self.m.coordinates)
        self.assertEqual(len(self.m.coordinates), 2)

    def test_query_edit(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        e = self.m.edits.edits[0]
        with mock.patch('builtins.input', return_value="y"):
            self.m.query_edit(e)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)

        with mock.patch('builtins.input', return_value="n"):
            self.m.query_edit(e)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)

    def test_query_edits(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        e = [self.m.edits.edits[0], self.m.edits.edits[1]]
        with mock.patch('builtins.input', return_value="y"):
            self.m.query_edits(e)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)

        with mock.patch('builtins.input', return_value="n"):
            self.m.query_edits(e)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)

    def test_add_query_to_edits(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        e = [self.m.edits.edits[0], self.m.edits.edits[2]]
        self.m.add_query_to_edits(e)
        self.assertEqual(self.m.edits.edits[0].edit_query, True)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, True)

    def test_find_common_edits_none(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.find_common_edits()
        self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[0].edit_query, False)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, False)

    def test_find_common_edits_last_2_common(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "c", "ref", 8))
        self.m.find_common_edits()
        assert self.m.edits.edits[0].edit_accepted == True
        assert self.m.edits.edits[1].edit_accepted == True
        assert self.m.edits.edits[2].edit_accepted == True
        assert self.m.edits.edits[3].edit_accepted == True
        assert self.m.edits.edits[0].edit_query == False
        assert self.m.edits.edits[1].edit_query == False
        assert self.m.edits.edits[2].edit_query == True
        assert self.m.edits.edits[3].edit_query == True

    def test_find_common_edits_last_2_common_interactive(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "c", "ref", 8))
        with mock.patch('builtins.input', return_value="y"):
            self.m.find_common_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)
        with mock.patch('builtins.input', return_value="n"):
            self.m.find_common_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_common_edits_last_2_common_min_3(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "c", "ref", 8))
        self.m.find_common_edits(min_occurrence=3)
        self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[0].edit_query, False)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, False)
        self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_similar_edits_none(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.find_similar_edits()
        self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[0].edit_query, False)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, False)

    def test_find_similar_edits_last_2_similar(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        self.m.find_similar_edits()
        self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[0].edit_query, False)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, True)
        self.assertEqual(self.m.edits.edits[3].edit_query, True)

    def test_find_similar_edits_last_2_similar_interactive_yy(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        mymock = mock.Mock()
        mymock.side_effect = ['y', 'y']
        with mock.patch('builtins.input', mymock):
            self.m.find_similar_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_similar_edits_last_2_similar_interactive_yn(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        mymock = mock.Mock()
        mymock.side_effect = ['y', 'n']
        with mock.patch('builtins.input', mymock):
            self.m.find_similar_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_similar_edits_last_2_similar_interactive_ny(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        mymock = mock.Mock()
        mymock.side_effect = ['n', 'y']
        with mock.patch('builtins.input', mymock):
            self.m.find_similar_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, False)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_similar_edits_last_2_similar_interactive_nn(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        mymock = mock.Mock()
        mymock.side_effect = ['n', 'n']
        with mock.patch('builtins.input', mymock):
            self.m.find_similar_edits(interactive=True)
            self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
            self.assertEqual(self.m.edits.edits[0].edit_query, False)
            self.assertEqual(self.m.edits.edits[1].edit_query, False)
            self.assertEqual(self.m.edits.edits[2].edit_query, False)
            self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_find_similar_edits_last_2_similar_min_3(self):
        self.m.edits = EditFile()
        self.m.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        self.m.edits.add_edit(Edit(0, 5, "g", "c", "ref", 6))
        self.m.edits.add_edit(Edit(0, 8, "g", "c", "ref", 8))
        self.m.edits.add_edit(Edit(1, 9, "g", "t", "ref", 8))
        self.m.find_similar_edits(min_occurrence=3)
        self.assertEqual(self.m.edits.edits[0].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[1].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[2].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[3].edit_accepted, True)
        self.assertEqual(self.m.edits.edits[0].edit_query, False)
        self.assertEqual(self.m.edits.edits[1].edit_query, False)
        self.assertEqual(self.m.edits.edits[2].edit_query, False)
        self.assertEqual(self.m.edits.edits[3].edit_query, False)

    def test_run_output_file_exists(self):
        self.assertRaises(AssertionError, self.m.run, data_dir, output_file="%s/consensus.fasta.edits" % data_dir)

    def test_run_output_file_defaults(self):
        tmp_file = os.path.join(data_dir, 'tmp.all.edits')
        self.m.run(data_dir, output_file="%s/all.edits" % data_dir)
        expect_file = os.path.join(data_dir, 'expected_edits')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)

    def test_run_output_file_tmp_exists(self):
        test_dir = "%s/tmp_exists" %data_dir
        tmp_file = os.path.join(test_dir, 'all.edits')
        self.m.run(test_dir, output_file=tmp_file)
        expect_file = os.path.join(test_dir, 'expected_edits')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)
        self.m.run(test_dir, output_file=tmp_file)

    def test_run_output_file_with_features(self):
        tmp_file = os.path.join(data_dir, 'tmp.all.edits')
        self.m.run(data_dir, output_file="%s/all.edits" % data_dir, features="ORF1")
        expect_file = os.path.join(data_dir, 'expected_edits_features')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)

    def test_run_output_file_save_output_file(self):
        tmp_file = os.path.join(data_dir, 'all.edits')
        self.m.run(data_dir, output_file="%s/all.edits" % data_dir, min_occurence=10)
        expect_file = os.path.join(data_dir, 'expected_edits_no_tmp_needed')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)
