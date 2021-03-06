import os
import unittest
import json
import filecmp

from genofunk.sequence_utils import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class TestSequenceUtils(unittest.TestCase):
    def test_get_coordinates_from_json_simple_pairs(self):
        json_value = {
          "start": 30,
          "end": 40,
          "strand": 1
        }
        coordinates = get_coordinates_from_json(json_value, pairs=True)
        expected = [[30, 40]]
        self.assertEqual(expected, coordinates)

    def test_get_coordinates_from_json_simple_no_pairs(self):
        json_value = {
          "start": 30,
          "end": 40,
          "strand": 1
        }
        coordinates = get_coordinates_from_json(json_value, pairs=False)
        expected = [30, 40]
        self.assertEqual(expected, coordinates)

    def test_get_coordinates_from_json_join_pairs(self):
        json_value = {
          "join": [
            { "start": 0, "end": 11, "strand": 1 },
            { "start": 10, "end": 20, "strand": 1 }
          ]
        }
        coordinates = get_coordinates_from_json(json_value, pairs=True)
        expected = [[0,11],[10,20]]
        self.assertEqual(expected, coordinates)

    def test_get_coordinates_from_json_join_no_pairs(self):
        json_value = {
          "join": [
            { "start": 0, "end": 11, "strand": 1 },
            { "start": 10, "end": 20, "strand": 1 }
          ]
        }
        coordinates = get_coordinates_from_json(json_value, pairs=False)
        expected = [0,11,10,20]
        self.assertEqual(expected, coordinates)

    def test_is_open_reading_frame_wrong_start(self):
        amino_acid_sequence = "NATIL*"
        result = is_open_reading_frame(amino_acid_sequence)
        self.assertFalse(result)

    def test_is_open_reading_frame_wrong_end(self):
        amino_acid_sequence = "MNATIL*S"
        result = is_open_reading_frame(amino_acid_sequence)
        self.assertFalse(result)

    def test_is_open_reading_frame_stop_in_middle(self):
        amino_acid_sequence = "MNATIL*S*"
        result = is_open_reading_frame(amino_acid_sequence, allow_stop_codons_in_middle=False)
        self.assertFalse(result)

    def test_is_open_reading_frame_stop_in_middle_allowed(self):
        amino_acid_sequence = "MNATIL*S*"
        result = is_open_reading_frame(amino_acid_sequence, allow_stop_codons_in_middle=True)
        self.assertTrue(result)