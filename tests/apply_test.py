import os
import unittest
from unittest import mock

import filecmp

from genofunk import apply
from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'apply')

class TestApply(unittest.TestCase):
    def setUp(self):
        self.a = apply.Apply()
        editfile = os.path.join(data_dir, 'consensus.fasta.edits')
        self.a.load_input_files(data_dir, editfile)

    def test_load_consensus_empty_directory(self):
        empty_dir = os.path.join(data_dir, "empty")
        self.assertRaises(AssertionError, self.a.load_consensus, empty_dir)

    def test_load_consensus_missing_pair(self):
        missing_dir = os.path.join(data_dir, "missing_consensus")
        self.assertRaises(AssertionError, self.a.load_consensus, missing_dir)

    def test_load_consensus(self):

        self.assertIsNotNone(self.a.consensus_sequence)
        self.assertEqual(len(self.a.consensus_sequence), 12)

    def test_load_input_files(self):
        self.assertIsNotNone(self.a.consensus_sequence)
        self.assertIsNotNone(self.a.edits)
        self.assertEqual(len(self.a.edits.edits), 14)
        self.assertEqual(len(self.a.consensus_sequence), 12)

    def test_apply_loaded_edits_default(self):
        pass

    def test_apply_loaded_edits_no_filter_by_accepted(self):
        pass

    def test_save_updated_consensuses(self):
        pass