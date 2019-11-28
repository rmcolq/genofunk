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
        if not os.path.exists(filepath):
            logging.error("Alignment file %s does not exist!" % filepath)
        self.amino_acid_alignment = AlignIO.read(open(filepath), filetype)
        logging.debug("The alignment file contains %d records" % len(self.amino_acid_alignment))
        assert (len(self.amino_acid_alignment) > 0)

    def load_input_files(self, consensus_filepath, alignment_filepath):
        self.load_consensus(consensus_filepath)
        self.load_alignment(alignment_filepath)
        logging.debug("Found %d consensus sequences and %d alignment sequences" %(int(len(self.consensus_sequences)), int(len(self.amino_acid_alignment))))
        if int(len(self.consensus_sequences)) != int(len(self.amino_acid_alignment)):
            logging.error("Alignment file and consensus file contain different numbers of sequences!")
            assert(len(self.consensus_sequences) == len(self.amino_acid_alignment))

    def generate_nucleotide_alignment(self):
        self.nucleotide_alignment = {}
        for record in self.amino_acid_alignment:
            nucleotide_sequence = self.consensus_sequences[record.id].seq
            logging.debug("3 x alignment length %d is %d, at least nucleotide sequence length %d" % (len(record.seq),
                                                                                            3*(len(record.seq)),
                                                                                            len(nucleotide_sequence)))
            no_dash = str(record.seq).replace("-","")
            logging.debug("3 x alignment length %d is %d, at least nucleotide sequence length %d" % (len(no_dash),
                                                                                                     3 * (len(no_dash)),
                                                                                              len(nucleotide_sequence)))
            assert 3*(len(record.seq)) >= len(nucleotide_sequence)
            aligned_nucleotide_sequence = ""
            pos = 0
            for i,letter in enumerate(record.seq):
                if letter == "-" and nucleotide_sequence[pos:pos+3] == "NNN":
                    aligned_nucleotide_sequence += "---"
                    pos += 3
                elif letter == "-":
                    aligned_nucleotide_sequence += "---"
                elif letter not in "*ACDEFGHIKLMNPQRSTVWXY":
                    print(letter)
                else:
                    aligned_nucleotide_sequence += nucleotide_sequence[pos:pos+3]
                    pos += 3
                #print(record.seq[:i+1])
                #print(aligned_nucleotide_sequence)
                #if i == 60:
                #    break
            self.nucleotide_alignment[record.id] = SeqRecord(aligned_nucleotide_sequence, id=record.id, description="")
            logging.debug("Have included %d bases of nucleotide sequnce of length %d for record %s" %(pos,
                                                                                  len(nucleotide_sequence), record.id))
            assert pos >= len(nucleotide_sequence)
            #break

    def save_nucleotide_alignment(self, filepath):
        SeqIO.write([self.nucleotide_alignment[i] for i in self.nucleotide_alignment], filepath, "fasta")

    def run(self, consensus_filepath, alignment_filepath, output_filepath):
        self.load_input_files(consensus_filepath, alignment_filepath)
        self.generate_nucleotide_alignment()
        self.save_nucleotide_alignment(output_filepath)
