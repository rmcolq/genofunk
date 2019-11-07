import os
import sys
import logging
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

class Translate():
    def __init__(self):
        self.consensus_sequence = None
        self.amino_acid_alignment = None
        self.nucleotide_alignment = None

    def load_consensus(self, filepath, filetype="fasta"):
        """
        Loads the (corrected nucleotide) consensus file
        :param filepath:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
        self.consensus_sequence = []

        logging.debug("Loading consensus %s %s" % (filetype, filepath))
        if not os.path.exists(filepath):
            logging.error("Consensus filepath %s does not exist!" % filepath)
        records = list(SeqIO.parse(filepath, filetype))
        logging.debug("The consensus file contains %d records" % len(records))
        assert (len(records) > 0)

    def load_alignment(self, filepath, filetype="fasta"):
        if not os.path.exists(filepath):
            logging.error("Alignment file %s does not exist!" % filepath)
        self.amino_acid_alignment = AlignIO.read(open(filepath), filetype)

    def load_input_files(self, consensus_filepath, alignment_filepath):
        self.load_consensus(consensus_filepath)
        self.load_alignment(alignment_filepath)
        if len(self.consensus_sequence) != len(self.amino_acid_alignment):
            logging.error("Alignment file and consensus file contain different numbers of sequences!")
            assert(len(self.consensus_sequence) == len(self.amino_acid_alignment))

    def generate_nucleotide_alignment(self):
        for record in self.amino_acid_alignment:
            print(record.id)
            print(self.consensus_sequence[record.id])


    def run(self, consensus_filepath, alignment_filepath):
        self.load_input_files(consensus_filepath, alignment_filepath)
        self.generate_nucleotide_alignment()
