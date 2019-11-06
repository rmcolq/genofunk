import os
import unittest
from Bio import SeqIO
import filecmp

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'editfile')

class TestEditFile(unittest.TestCase):

    def setUp(self):
        consensus_filepath = os.path.join(data_dir, 'consensus.fasta')
        self.consensus_sequence = list(SeqIO.parse(consensus_filepath, 'fasta'))

    def test_edit_attributes(self):
        sequence_id = 8
        sequence_position = 3
        edit_from = ""
        edit_to = "N"
        reference_id = "ref"
        reference_position = 5
        e = Edit(sequence_id, sequence_position, edit_from, edit_to, reference_id, reference_position)
        self.assertEqual(e.sequence_id, sequence_id)
        self.assertEqual(e.sequence_position, sequence_position)
        self.assertEqual(e.edit_from, edit_from)
        self.assertEqual(e.edit_to, edit_to)
        self.assertEqual(e.reference_id, reference_id)
        self.assertEqual(e.reference_position, reference_position)
        self.assertEqual(e.edit_applied, False)
        self.assertEqual(e.edit_accepted, True)


    def test_edit_attributes_with_accepted(self):
        sequence_id = 8
        sequence_position = 3
        edit_from = ""
        edit_to = "N"
        reference_id = "ref"
        reference_position = 5
        edit_accepted = False
        e = Edit(sequence_id, sequence_position, edit_from, edit_to, reference_id, reference_position, edit_accepted)
        self.assertEqual(e.sequence_id, sequence_id)
        self.assertEqual(e.sequence_position, sequence_position)
        self.assertEqual(e.edit_from, edit_from)
        self.assertEqual(e.edit_to, edit_to)
        self.assertEqual(e.reference_id, reference_id)
        self.assertEqual(e.reference_position, reference_position)
        self.assertEqual(e.edit_applied, False)
        self.assertEqual(e.edit_accepted, False)

    def test_edit_missing_attribute(self):
        sequence_id = 8
        sequence_position = None
        edit_from = ""
        edit_to = "N"
        reference_id = "ref"
        reference_position = 5
        self.assertRaises(TypeError, Edit, sequence_id, edit_from, edit_to, reference_id, reference_position)

    def test_edit_equals(self):
        e1 = Edit(0, 0, "a", "g", "ref", 1)
        e2 = Edit(1, 0, "a", "g", "ref", 1)
        e3 = Edit(0, 1, "a", "g", "ref", 1)
        e4 = Edit(0, 0, "c", "g", "ref", 1)
        e5 = Edit(0, 0, "a", "t", "ref", 1)
        e6 = Edit(0, 0, "a", "g", "reference", 1)
        e7 = Edit(0, 0, "a", "g", "ref", 2)
        e8 = Edit(0, 0, "a", "g", "ref", 1, edit_accepted=False)
        e9 = Edit(0, 0, "a", "g", "ref", 1, edit_applied=True)

        self.assertEqual(e1, e1)
        self.assertEqual(e2, e2)
        self.assertEqual(e3, e3)
        self.assertEqual(e4, e4)
        self.assertEqual(e5, e5)
        self.assertEqual(e6, e6)
        self.assertEqual(e7, e7)
        self.assertEqual(e8, e8)
        self.assertEqual(e9, e9)
        self.assertEqual(e1, e8)
        self.assertEqual(e8, e1)
        self.assertEqual(e1, e9)
        self.assertEqual(e9, e1)
        self.assertEqual(e8, e9)
        self.assertEqual(e9, e8)

    def test_edit_not_equals(self):
        e1 = Edit(0, 0, "a", "g", "ref", 1)
        e2 = Edit(1, 0, "a", "g", "ref", 1)
        e3 = Edit(0, 1, "a", "g", "ref", 1)
        e4 = Edit(0, 0, "c", "g", "ref", 1)
        e5 = Edit(0, 0, "a", "t", "ref", 1)
        e6 = Edit(0, 0, "a", "g", "reference", 1)
        e7 = Edit(0, 0, "a", "g", "ref", 2)
        e8 = Edit(0, 0, "a", "g", "ref", 1, edit_accepted=False)
        e9 = Edit(0, 0, "a", "g", "ref", 1, edit_applied=True)

        self.assertNotEqual(e1, e2)
        self.assertNotEqual(e1, e3)
        self.assertNotEqual(e1, e4)
        self.assertNotEqual(e1, e5)
        self.assertNotEqual(e1, e6)
        self.assertNotEqual(e1, e7)
        self.assertNotEqual(e2, e3)
        self.assertNotEqual(e2, e4)
        self.assertNotEqual(e2, e5)
        self.assertNotEqual(e2, e6)
        self.assertNotEqual(e2, e7)
        self.assertNotEqual(e2, e8)
        self.assertNotEqual(e2, e9)
        self.assertNotEqual(e3, e4)
        self.assertNotEqual(e3, e5)
        self.assertNotEqual(e3, e6)
        self.assertNotEqual(e3, e7)
        self.assertNotEqual(e3, e8)
        self.assertNotEqual(e3, e9)
        self.assertNotEqual(e4, e5)
        self.assertNotEqual(e4, e6)
        self.assertNotEqual(e4, e7)
        self.assertNotEqual(e4, e8)
        self.assertNotEqual(e4, e9)
        self.assertNotEqual(e5, e6)
        self.assertNotEqual(e5, e7)
        self.assertNotEqual(e5, e8)
        self.assertNotEqual(e5, e9)
        self.assertNotEqual(e6, e7)
        self.assertNotEqual(e6, e8)
        self.assertNotEqual(e6, e9)
        self.assertNotEqual(e7, e8)
        self.assertNotEqual(e7, e9)

    def test_edit_less_than(self):
        e1 = Edit(0, 0, "a", "g", "ref", 1)
        e2 = Edit(1, 0, "a", "g", "ref", 1)
        e3 = Edit(0, 1, "a", "g", "ref", 1)
        e4 = Edit(0, 0, "c", "g", "ref", 1)
        e5 = Edit(0, 0, "a", "t", "ref", 1)
        e6 = Edit(0, 0, "a", "g", "reference", 1)
        e7 = Edit(0, 0, "a", "g", "ref", 2)
        e8 = Edit(0, 0, "a", "g", "ref", 1, edit_accepted=False)
        e9 = Edit(0, 0, "a", "g", "ref", 1, edit_applied=True)

        self.assertTrue(e1 < e2)
        self.assertTrue(e1 < e3)
        self.assertTrue(e1 < e4)
        self.assertTrue(e1 < e5)
        self.assertTrue(e1 < e6)
        self.assertTrue(e1 < e7)
        self.assertFalse(e1 < e8)
        self.assertFalse(e1 < e9)

    def test_edit_apply_insertion(self):
        record_id = 0
        e = Edit(record_id, 3, '', "N", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], 0)
        self.assertEqual(str(self.consensus_sequence[record_id].seq),
                         "attNaacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(e.edit_from, "")
        self.assertEqual(e.edit_to, "N")

    def test_edit_apply_deletion(self):
        record_id = 0
        e = Edit(record_id, 3, 'N', "", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], 0)
        self.assertEqual(str(self.consensus_sequence[record_id].seq),
                         "attacgcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(e.edit_from, "a")
        self.assertEqual(e.edit_to, "")

    def test_edit_apply_deletion_with_offset(self):
        record_id = 0
        e = Edit(record_id, 3, 'N', "", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], 2)
        self.assertEqual(str(self.consensus_sequence[record_id].seq),
                         "attaagcgcatctggaaattaacacccatgaaggccgcaacgatacccatgaacgcgaactgattgtggaagatgcgcatattacctaa")
        self.assertEqual(e.edit_from, "c")
        self.assertEqual(e.edit_to, "")

    def test_edit_invertible(self):
        record_id = 0
        original = self.consensus_sequence[record_id]

        e = Edit(record_id, 3, 'N', "", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], 0)
        e.remove_edit(self.consensus_sequence[record_id], 0)
        self.assertEqual(str(original.seq), str(self.consensus_sequence[record_id].seq))

        e = Edit(record_id, 3, '', "N", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], 1)
        e.remove_edit(self.consensus_sequence[record_id], 1)
        self.assertEqual(str(original.seq), str(self.consensus_sequence[record_id].seq))

        e = Edit(record_id, 3, '', "N", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id], -1)
        e.remove_edit(self.consensus_sequence[record_id], -1)
        self.assertEqual(str(original.seq), str(self.consensus_sequence[record_id].seq))

    def test_editfile_equals(self):
        file1 = os.path.join(data_dir, 'editfile1.csv')
        file2 = os.path.join(data_dir, 'editfile2.csv')

        editfile1 = EditFile(file1)
        editfile2 = EditFile(file2)

        self.assertEqual(editfile1, editfile1)
        self.assertEqual(editfile2, editfile2)
        self.assertNotEqual(editfile1, editfile2)
        self.assertNotEqual(editfile2, editfile1)

    def test_editfile_append_same_as_load(self):
        editfile = EditFile()
        file1 = os.path.join(data_dir, 'editfile1.csv')
        editfile.append(file1)

        expect_editfile = EditFile(file1)
        self.assertEqual(editfile, expect_editfile)

    def test_editfile_append_two_files(self):
        file1 = os.path.join(data_dir, 'editfile1.csv')
        file2 = os.path.join(data_dir, 'editfile2.csv')
        editfile = EditFile(file1)
        editfile.append(file2)

        expect_file = os.path.join(data_dir, 'editfile1and2.csv')
        expect_editfile = EditFile(expect_file)
        self.assertEqual(editfile, expect_editfile)

    def test_editfile_sort(self):
        file = os.path.join(data_dir, 'editfile1and2.csv')
        editfile = EditFile(file)
        editfile.sort()

        expect_file = os.path.join(data_dir, 'editfile_sorted.csv')
        expect_editfile = EditFile(expect_file)
        self.assertEqual(editfile, expect_editfile)

    def test_editfile_save(self):
        editfile = EditFile()
        record_id = 0
        e = Edit(record_id, 3, 'N', "", "hobbit", 5)
        e.apply_edit(self.consensus_sequence[record_id])
        editfile.add_edit(e)
        e = Edit(record_id, 6, 'N', "", "hobbit", 10)
        editfile.add_edit(e)


        tmp_file = os.path.join(data_dir, 'tmp_saved.csv')
        expect_file = os.path.join(data_dir, 'expect_saved.csv')
        editfile.save(tmp_file)
        self.assertTrue(filecmp.cmp(tmp_file, expect_file, shallow=False))
        os.unlink(tmp_file)