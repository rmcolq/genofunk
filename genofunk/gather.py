import os
import logging
import glob
from Bio import SeqIO

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

class Gather():
    def __init__(self):
        self.consensus_sequence = None
        self.edits = None

    def load_from_directory(self, directory, filetype="fasta"):
        self.consensus_sequence = []
        self.edits = EditFile()

        edit_files = glob.iglob("%s/*.edits" %directory)
        for edit_file in edit_files:
            consensus_file = edit_file.replace(".edits","")
            if not os.path.exists(consensus_file):
                logging.error("Paired consensus file %s does not exist!" % consensus_file)

            self.consensus_sequence.extend(list(SeqIO.parse(consensus_file, filetype)))
            self.edits.append(edit_file)
            logging.debug("Now have %d consensus records and %d edits" %(len(self.consensus_sequence), len(self.edits.edits)))

    def run(self, directory):
        self.load_from_directory(directory)