import os
import unittest
from unittest import mock

import filecmp

from genofunk import gather
from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'gather')

class TestGather(unittest.TestCase):
    def setUp(self):
        self.g = gather.Gather()

    def test_load_from_directory_empty(self):
        empty_dir = os.path.join(data_dir, "empty")
        self.assertRaises(AssertionError, self.g.load_from_directory, empty_dir)

    def test_load_from_directory_missing_pair(self):
        missing_dir = os.path.join(data_dir, "missing_consensus")
        self.assertRaises(AssertionError, self.g.load_from_directory, missing_dir)

    def test_load_from_directory(self):
        self.g.load_from_directory(data_dir)
        self.assertIsNotNone(self.g.edits)
        self.assertIsNotNone(self.g.consensus_sequence)
        self.assertEqual(len(self.g.edits.edits), 14)
        self.assertEqual(len(self.g.consensus_sequence), 12)

    def test_query_edit(self):
        self.g.edits = EditFile()
        self.g.edits.add_edit(Edit(0, 0, "a", "g", "ref", 1))
        e = self.g.edits.edits[0]
        print(e)
        with mock.patch('builtins.input', return_value="y"):
            self.g.query_edit(e)
            assert e.edit_accepted == False

        with mock.patch('builtins.input', return_value="n"):
            self.g.query_edit(e)
            assert e.edit_accepted == True