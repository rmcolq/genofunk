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

    def test_apply_loaded_edits_edits_none(self):
        a = apply.Apply()
        a.load_consensus(data_dir)
        a.apply_loaded_edits()
        for name in a.consensus_sequence:
            self.assertEqual(a.consensus_sequence[name].seq, self.a.consensus_sequence[name].seq)

    def test_apply_loaded_edits_edits_empty(self):
        a = apply.Apply()
        a.load_consensus(data_dir)
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
                         "attaacgcgcatctggaaattaacacccatgaaggcccaacgatacccatgaacgcgaactgatttgggaagatgcgcatattacctaa")
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
                         "attaacgcgcatctggattaacacccatgaaggccccaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(str(self.a.consensus_sequence["seq2_with_3_insertions"].seq),
                         "aacaccgcgaacgcgagcacctatgatattcgcacctattgggaaacccatctggaatttattctgctggaagattggattacccatacccatgaagaaaacgatagcttttggcgcatgagctaa")
        num_applied = 0
        for edit in self.a.edits.edits:
            if edit.edit_applied:
                num_applied += 1
        self.assertEqual(num_applied, 12)

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

    def test_save_updated_consensuses_na(self):
        self.a.apply_loaded_edits(filter_by_accepted=False)
        tmp_file = os.path.join(data_dir, 'updated_aa_consensus.fasta')
        self.a.save_updated_consensuses(filepath=tmp_file, amino_acid=True)
        expect_file = os.path.join(data_dir, 'expect_updated_aa_consensus.fasta')
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)
