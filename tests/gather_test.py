import os
import unittest
import filecmp

from genofunk import gather

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