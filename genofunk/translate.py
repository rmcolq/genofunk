import os
import sys
import logging
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Translate():
    def __init__(self):
        self.consensus_sequences = None
        self.amino_acid_alignment = None
        self.nucleotide_alignment = None

    def load_consensus(self, filepath, filetype="fasta"):
        """
        Loads the (corrected nucleotide) consensus file
        :param filepath:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
        self.consensus_sequences = []

        logging.debug("Loading consensus %s %s" % (filetype, filepath))
        if not os.path.exists(filepath):
            logging.error("Consensus filepath %s does not exist!" % filepath)
        self.consensus_sequences = SeqIO.index(filepath, filetype)
        logging.debug("The consensus file contains %d records" % len(self.consensus_sequences))
        assert (len(self.consensus_sequences) > 0)

    def load_alignment(self, filepath, filetype="fasta"):
        """
        Loads the amino acid alignment
        :param filepath:
        :param filetype: (Default: FASTA)
        :return:
        """
        if not os.path.exists(filepath):
            logging.error("Alignment file %s does not exist!" % filepath)
        self.amino_acid_alignment = AlignIO.read(open(filepath), filetype)
        logging.debug("The alignment file contains %d records" % len(self.amino_acid_alignment))
        assert (len(self.amino_acid_alignment) > 0)

    def load_input_files(self, consensus_filepath, alignment_filepath):
        """
        Loads both nucleotide consensus file and amino acid multiple sequence alignment file and checks have same number
        of sequences
        :param consensus_filepath:
        :param alignment_filepath:
        :return:
        """
        self.load_consensus(consensus_filepath)
        self.load_alignment(alignment_filepath)
        logging.debug("Found %d consensus sequences and %d alignment sequences" %(int(len(self.consensus_sequences)),
                                                                                  int(len(self.amino_acid_alignment))))
        if int(len(self.consensus_sequences)) != int(len(self.amino_acid_alignment)):
            logging.error("Alignment file and consensus file contain different numbers of sequences!")
            assert(len(self.consensus_sequences) == len(self.amino_acid_alignment))

    def generate_nucleotide_alignment(self):
        """
        Uses amino acid alignment to guide a nucleotide alignment
        :return:
        """
        self.nucleotide_alignment = {}
        for record in self.amino_acid_alignment:
            nucleotide_sequence = self.consensus_sequences[record.id].seq
            assert 3*(len(record.seq)) >= len(nucleotide_sequence)
            aligned_nucleotide_sequence = ""
            pos = 0
            for i,letter in enumerate(record.seq):
                if letter == "-":
                    aligned_nucleotide_sequence += "---"
                elif letter not in "*ACDEFGHIKLMNPQRSTVWXY":
                    print(letter)
                else:
                    while pos+3 > len(nucleotide_sequence):
                        nucleotide_sequence += "-"
                    aligned_nucleotide_sequence += nucleotide_sequence[pos:pos+3]
                    pos += 3
            self.nucleotide_alignment[record.id] = SeqRecord(aligned_nucleotide_sequence, id=record.id, description="")
            logging.debug("Have included %d bases of nucleotide sequnce of length %d for record %s" %(pos,
                                                                                  len(nucleotide_sequence), record.id))
            assert pos >= len(nucleotide_sequence)
            #break

    def check_nucleotide_alignment_lengths(self):
        """
        Checks all inferred nucleotide alignment sequences have the same length
        :return:
        """
        lengths = list(set([len(self.nucleotide_alignment[i]) for i in self.nucleotide_alignment]))
        if len(lengths) > 1:
            logging.debug("Aligned nucleotide sequence lengths: %s" %",".join(lengths))
            assert len(lengths) == 1

    def save_nucleotide_alignment(self, filepath):
        """
        Saves nucleotide alignment
        :param filepath:
        :return:
        """
        SeqIO.write([self.nucleotide_alignment[i] for i in self.nucleotide_alignment], filepath, "fasta")

    def run(self, consensus_filepath, alignment_filepath, output_filepath):
        """
        Uses an amino acid multiple sequence alignment to guide creation of a nucleotide multiple sequence alignment
        :param consensus_filepath:
        :param alignment_filepath:
        :param output_filepath:
        :return:
        """
        self.load_input_files(consensus_filepath, alignment_filepath)
        self.generate_nucleotide_alignment()
        self.check_nucleotide_alignment_lengths()
        self.save_nucleotide_alignment(output_filepath)
        self.consensus_sequences.close()
