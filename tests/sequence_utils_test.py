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
        coordinates = get_coordinates_from_json(json_value, pairs=True)
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