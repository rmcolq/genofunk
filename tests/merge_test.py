import os
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