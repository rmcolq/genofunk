import os
import unittest
import filecmp
import json

from genofunk import apply
from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'apply')

class TestApply(unittest.TestCase):
    def setUp(self):
        self.a = apply.Apply()
        editfile = os.path.join(data_dir, 'all.edits')
        self.a.load_input_files(data_dir, editfile)

    def test_load_coordinates_file_no_coordinates(self):
        a = apply.Apply()
        coordinate_file = os.path.join(data_dir, "doesnt_exist.coordinates")
        self.assertRaises(AssertionError, a.load_coordinates_file, coordinate_file)

    def test_load_coordinates_file_does_not_exist(self):
        a = apply.Apply()
        a.coordinates = {}
        coordinate_file = os.path.join(data_dir, "doesnt_exist.coordinates")
        self.assertRaises(AssertionError, a.load_coordinates_file, coordinate_file)

    def test_load_coordinates_file_empty(self):
        a = apply.Apply()
        a.coordinates = {}
        coordinate_file = os.path.join(data_dir, "empty.coordinates")
        self.assertRaises(json.JSONDecodeError, a.load_coordinates_file, coordinate_file)

    def test_load_coordinates_file_no_feature_list(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = None
        a = apply.Apply()
        a.coordinates = {}
        a.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "ORF2": {"seq2": {"start": 0, "end": 126}, "seq3": {"start": 10, "end": 136}}}
        self.assertEqual(a.coordinates, expect_coordinates)

    def test_load_coordinates_file_with_feature_list(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = ["ORF1", "idontexistinasample"]
        a = apply.Apply()
        a.coordinates = {}
        a.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        self.assertEqual(a.coordinates, expect_coordinates)

    def test_load_coordinates_file_two_files(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        features_list = ["ORF1", "ORF2", "idontexistinasample"]
        a = apply.Apply()
        a.coordinates = {}
        a.load_coordinates_file(coordinate_file, features_list)
        coordinate_file = os.path.join(data_dir, "simple2.coordinates")
        a.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100},
                                       "seq4": {"start": 0, "end": 126}, "seq5": {"start": 10, "end": 136}},
                              "ORF2": {"seq2": {"start": 0, "end": 126}, "seq3": {"start": 10, "end": 136}},
                              "idontexistinasample": {}}
        self.assertEqual(a.coordinates, expect_coordinates)

    def test_load_coordinates_file_with_join(self):
        coordinate_file = os.path.join(data_dir, "simple2.coordinates")
        features_list = ["ORF3", "idontexistinasample"]
        a = apply.Apply()
        a.coordinates = {}
        a.load_coordinates_file(coordinate_file, features_list)
        expect_coordinates = {"ORF3": {"seq6": {"join": [{"start": 0, "end": 124}, {"start": 125, "end": 128}]}},
                              "idontexistinasample": {}}
        self.assertEqual(a.coordinates, expect_coordinates)

    def test_load_coordinates_file_repeat_ids(self):
        coordinate_file = os.path.join(data_dir, "simple.coordinates")
        self.assertRaises(AssertionError, self.a.load_coordinates_file, coordinate_file)

    def test_load_consensus_file_no_consensus(self):
        a = apply.Apply()
        consensus_file = os.path.join(data_dir, "doesnt_exist.fasta")
        self.assertRaises(AssertionError, a.load_consensus_file, consensus_file)

    def test_load_consensus_file_does_not_exist(self):
        a = apply.Apply()
        a.consensus = {}
        consensus_file = os.path.join(data_dir, "doesnt_exist.fasta")
        self.assertRaises(AssertionError, a.load_consensus_file, consensus_file)

    def test_load_consensus_file_empty(self):
        a = apply.Apply()
        a.consensus = {}
        consensus_file = os.path.join(data_dir, "empty.fasta")
        self.assertRaises(AssertionError, a.load_consensus_file, consensus_file)

    def load_consensus_file(self):
        a = apply.Apply()
        a.consensus = {}
        consensus_file = os.path.join(data_dir, "consensus.fasta")
        a.load_consensus_file(consensus_file)
        self.assertEqual(len(a.consensus_sequence), 12)

    def load_consensus_file_repeat_ids(self):
        a = apply.Apply()
        a.consensus = {}
        consensus_file = os.path.join(data_dir, "consensus.fasta")
        a.load_consensus_file(consensus_file)
        self.assertRaises(AssertionError, a.load_consensus_file, consensus_file)


    def test_load_edits_in_range_no_coordinates(self):
        a = apply.Apply()
        edit_file = os.path.join(data_dir, "simple_edits")
        self.assertRaises(AssertionError, a.load_edits_in_range, edit_file)

    def test_load_edits_in_range_no_feature_list(self):
        a = apply.Apply()
        edit_file = os.path.join(data_dir, "simple_edits")
        features_list = None
        a.edits = EditFile()
        a.coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        a.load_edits_in_range(edit_file, features_list)
        expect_edits = EditFile(edit_file)
        self.assertEqual(a.edits, expect_edits)

    def test_load_edits_in_range_with_feature_list(self):
        a = apply.Apply()
        edit_file = os.path.join(data_dir, "simple_edits")
        features_list = ["ORF1", "idontexistinasample"]
        a.coordinates = {"ORF1": {"seq1": {"start": 0, "end": 90}, "seq3": {"start": 10, "end": 100}},
                              "idontexistinasample": {}}
        a.load_edits_in_range(edit_file, features_list)
        expect_edit_file = os.path.join(data_dir, "simple_edits_expect")
        expect_edits = EditFile(expect_edit_file)
        self.assertEqual(a.edits, expect_edits)

    def test_load_input_files_empty(self):
        a = apply.Apply()
        edit_file = os.path.join(data_dir, "all.edits")
        empty_dir = os.path.join(data_dir, "empty")
        self.assertRaises(AssertionError, a.load_input_files, empty_dir, edit_file)

    def test_load_input_files_missing_pair(self):
        a = apply.Apply()
        edit_file = os.path.join(data_dir, "all.edits")
        missing_dir = os.path.join(data_dir, "missing_consensus")
        self.assertRaises(AssertionError, a.load_input_files, missing_dir, edit_file)

    def test_load_input_files(self):
        self.assertIsNotNone(self.a.edits)
        self.assertIsNotNone(self.a.coordinates)
        self.assertIsNotNone(self.a.consensus_sequence)
        self.assertEqual(len(self.a.edits.edits), 14)
        self.assertEqual(len(self.a.coordinates), 2)
        self.assertEqual(len(self.a.consensus_sequence), 12)

    def test_apply_loaded_edits_edits_empty(self):
        a = apply.Apply()
        editfile = os.path.join(data_dir, "all.edits")
        a.load_input_files(data_dir, edit_filepath=editfile)
        a.edits = EditFile()
        a.apply_loaded_edits()
        for name in a.consensus_sequence:
            self.assertEqual(a.consensus_sequence[name].seq, self.a.consensus_sequence[name].seq)

    def test_apply_loaded_edits_default(self):
        self.a.apply_loaded_edits()
        print(self.a.edits)
        self.assertEqual(str(self.a.consensus_sequence["seq1"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq2"].seq),
                         "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion_and_insertion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaactgatttgggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion_and_mismatch"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaacagattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_2_deletions"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaacNNattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_insertion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_2_insertions"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_mismatches"].seq),
                         "attaacgcgcatctggaaattatcacccatgaaggccgcaacgatacccatgatggcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_mismatch_insertion_deletion_mismatch"].seq),
                         "attaacgcgcatctggtaattaacacccatgaNggccgcaacgatacccatgaacgcgaagtgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_double_deletion_mismatch_insertion"].seq),
                         "attaacgcgcatctggNNattaacacccatgaaggccccaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq2_with_3_insertions"].seq),
                         "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        num_applied = 0
        for edit in self.a.edits.edits:
            if edit.edit_applied:
                num_applied += 1
        self.assertEqual(num_applied, 14)

    def test_apply_loaded_edits_no_filter_by_accepted(self):
        self.a.apply_loaded_edits(filter_by_accepted=False)
        print(self.a.edits)
        self.assertEqual(str(self.a.consensus_sequence["seq1"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq2"].seq),
                         "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion_and_insertion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaactgatttgggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_deletion_and_mismatch"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccNcaacgatacccatgaacgcgaacagattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_2_deletions"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaacNNattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_1_insertion"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_2_insertions"].seq),
                         "attaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_mismatches"].seq),
                         "attaacgcgcatctggaaattatcacccatgaaggccgcaacgatacccatgatggcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_mismatch_insertion_deletion_mismatch"].seq),
                         "attaacgcgcatctggtaattaacacccatgaNggccgcaacgatacccatgaacgcgaagtgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq1_with_double_deletion_mismatch_insertion"].seq),
                         "attaacgcgcatctggNNattaacacccatgaaggccccaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq2_with_3_insertions"].seq),
                         "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        for edit in self.a.edits.edits:
            self.assertEqual(edit.edit_applied, True)

    def test_save_updated_consensuses_na(self):
        self.a.apply_loaded_edits(filter_by_accepted=False)
        tmp_file = os.path.join(data_dir, 'updated_na_consensus.fasta')
        self.a.save_updated_consensuses(filepath=tmp_file)
        expect_file = os.path.join(data_dir, 'expect_updated_na_consensus.fasta')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)

    def test_save_updated_consensuses_aa(self):
        self.a.apply_loaded_edits(filter_by_accepted=False)
        tmp_file = os.path.join(data_dir, 'updated_aa_consensus.fasta')
        self.a.save_updated_consensuses(filepath=tmp_file, amino_acid=True)
        expect_file = os.path.join(data_dir, 'expect_updated_aa_consensus.fasta')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)

    def test_run_no_features(self):
        a = apply.Apply()
        editfile = os.path.join(data_dir, 'all.edits')
        tmp_prefix = os.path.join(data_dir, 'tmp')
        a.run(data_dir, editfile, tmp_prefix, features="")
        tmp_na_file = os.path.join(data_dir, 'tmp.na.fasta')
        expect_na_file = os.path.join(data_dir, 'expect_updated_na_consensus.fasta')
        self.assertTrue(filecmp.cmp(tmp_na_file, expect_na_file, shallow=False))
        os.unlink(tmp_na_file)
        tmp_aa_file = os.path.join(data_dir, 'tmp.aa.fasta')
        expect_aa_file = os.path.join(data_dir, 'expect_updated_aa_consensus.fasta')
        self.assertTrue(filecmp.cmp(tmp_aa_file, expect_aa_file, shallow=False))
        os.unlink(tmp_aa_file)

    def test_run_with_features(self):
        a = apply.Apply()
        editfile = os.path.join(data_dir, 'all.edits')
        tmp_prefix = os.path.join(data_dir, 'tmp')
        a.run(data_dir, editfile, tmp_prefix, features="ORF1")
        tmp_na_file = os.path.join(data_dir, 'tmp.ORF1.na.fasta')
        expect_na_file = os.path.join(data_dir, 'expect_updated_na_consensus_ORF1.fasta')
        self.assertTrue(filecmp.cmp(tmp_na_file, expect_na_file, shallow=False))
        os.unlink(tmp_na_file)
        tmp_aa_file = os.path.join(data_dir, 'tmp.ORF1.aa.fasta')
        expect_aa_file = os.path.join(data_dir, 'expect_updated_aa_consensus_ORF1.fasta')
        self.assertTrue(filecmp.cmp(tmp_aa_file, expect_aa_file, shallow=False))
        os.unlink(tmp_aa_file)

    def test_run_with_features_concat(self):
        a = apply.Apply()
        editfile = os.path.join(data_dir, 'all.edits')
        tmp_prefix = os.path.join(data_dir, 'tmp')
        a.run(data_dir, editfile, tmp_prefix, concat=True, features="ORF1,ORF2")
        tmp_na_file = os.path.join(data_dir, 'tmp.na.fasta')
        expect_na_file = os.path.join(data_dir, 'expect_updated_na_consensus_concat.fasta')
        self.assertTrue(filecmp.cmp(tmp_na_file, expect_na_file, shallow=False))
        os.unlink(tmp_na_file)
        tmp_aa_file = os.path.join(data_dir, 'tmp.aa.fasta')
        expect_aa_file = os.path.join(data_dir, 'expect_updated_aa_consensus_concat.fasta')
        self.assertTrue(filecmp.cmp(tmp_aa_file, expect_aa_file, shallow=False))
        os.unlink(tmp_aa_file)

    def test_run_with_no_features_concat(self):
        a = apply.Apply()
        editfile = os.path.join(data_dir, 'all.edits')
        tmp_prefix = os.path.join(data_dir, 'tmp')
        a.run(data_dir, editfile, tmp_prefix, concat=True)
        tmp_na_file = os.path.join(data_dir, 'tmp.na.fasta')
        expect_na_file = os.path.join(data_dir, 'expect_updated_na_consensus_concat.fasta')
        self.assertTrue(filecmp.cmp(tmp_na_file, expect_na_file, shallow=False))
        os.unlink(tmp_na_file)
        tmp_aa_file = os.path.join(data_dir, 'tmp.aa.fasta')
        expect_aa_file = os.path.join(data_dir, 'expect_updated_aa_consensus_concat.fasta')
        self.assertTrue(filecmp.cmp(tmp_aa_file, expect_aa_file, shallow=False))
        os.unlink(tmp_aa_file)
