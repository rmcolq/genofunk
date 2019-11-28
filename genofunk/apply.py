import os
import sys
import logging
import glob
from Bio import SeqIO
from Bio.Seq import Seq

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

class Apply():
    def __init__(self):
        self.consensus_sequence = None
        self.edits = None

    def load_consensus(self, directory, edit_filepath=None, filetype="fasta"):
        """
        Looks for pairs of *.fasta, *.fasta.edit files in a directory, and loads the *.fasta them into a single big
         list of consensus sequences
        :param directory:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
        self.consensus_sequence = {}

        edit_files = glob.glob("%s/*.edits" %directory)
        if edit_filepath:
            edit_files = list(filter(lambda x: not x.endswith(edit_filepath), edit_files))
        if len(edit_files) == 0:
            logging.error("No edit files found in directory %s" %directory)
            assert(len(edit_files) > 0)
        for edit_file in edit_files:
            consensus_file = edit_file.replace(".edits","")
            if not os.path.exists(consensus_file):
                logging.error("Paired consensus file %s does not exist!" % consensus_file)
                assert(os.path.exists(consensus_file))
            logging.debug("Loading consensus file %s and edit file %s" %(consensus_file, edit_file))
            record_dict = SeqIO.index(consensus_file, filetype)
            self.consensus_sequence.update(record_dict)
            logging.debug("Now have %d consensus records" %len(self.consensus_sequence))


    def load_input_files(self, directory, edit_filepath):
        self.load_consensus(directory, edit_filepath=edit_filepath)
        self.edits = EditFile(edit_filepath)
        self.edits.sort(reverse=True)
        logging.debug("Now have %d sorted edits" % len(self.edits.edits))

    def apply_loaded_edits(self, filter_by_accepted=True):
        """
        Apply the edits to the consensus nucleotide sequences in place (in reverse order to avoid offset errors)
        :return:
        """
        if not self.edits or len(self.edits.edits) == 0:
            return
        for edit in self.edits.edits:
            record = self.consensus_sequence[edit.sequence_id]
            edit.apply_edit(record, filter_by_accepted=filter_by_accepted)

    def save_updated_consensuses(self, filepath=None, coordinates=None, amino_acid=False):
        if filepath:
            out_handle = open(filepath, 'w')
        else:
            out_handle = sys.stdout

        for seq_name in self.consensus_sequence:
            seq = self.consensus_sequence[seq_name]
            assert (type(seq), Seq)
            if coordinates:
                (start, end) = coordinates
                seq = seq[start:end]
            if amino_acid:
                if len(seq) % 3 == 1:
                    seq = seq + "NN"
                elif len(seq) % 3 == 2:
                    seq = seq + "N"
                new_seq = seq.translate(stop_symbol='X')
                new_seq.id = seq.id
                new_seq.name = seq.name
                new_seq.description = seq.description
                seq = new_seq
                assert (len(seq)-1)*3 <= len(self.consensus_sequence[seq_name]) <= len(seq)*3
            SeqIO.write(seq, out_handle, "fasta")

        if filepath:
            out_handle.close()

    def run(self, directory, edit_filepath, output_prefix):
        self.load_input_files(directory, edit_filepath)
        self.apply_loaded_edits()
        self.save_updated_consensuses("%s.na.fasta" % output_prefix)
        self.save_updated_consensuses("%s.aa.fasta" % output_prefix, amino_acid=True)
