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

    def test_edit_missing_attribute(self):
        sequence_id = 8
        sequence_position = None
        edit_from = ""
        edit_to = "N"
        reference_id = "ref"
        reference_position = 5
        self.assertRaises(TypeError, Edit, sequence_id, edit_from, edit_to, reference_id, reference_position)

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