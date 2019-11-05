import os
import logging
from Bio import SeqIO
import json
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import parasail

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit


class Gather():
    def __init__(self, directory):
        self.consensus_sequence = None
        self.edits = None

    def load_reference_info(self, filepath):
        """
        Load the reference JSON, checking it contains the appropriate fields
        :param filepath:
        :return: dict (from JSON)
        """
        logging.debug("Loading reference JSON %s" % filepath)
        if not os.path.exists(filepath):
            logging.error("Reference filepath %s does not exist!" % filepath)
        with open(filepath) as json_file:
            data = json.load(json_file)
        logging.debug("Checking that JSON has correct format and contains the appropriate fields (sequence and orf) "
                      "for accession %s" % self.closest_accession)
        assert ('references' in data.keys())
        assert (self.closest_accession in data['references'].keys())
        assert ('sequence' in data['references'][self.closest_accession].keys())
        assert ('orf' in data['references'][self.closest_accession].keys())
        return data

    def load_consensus_sequence(self, filepath, filetype="fasta"):
        """
        Load consensus FASTA (or other file type) and check there is at least one record
        :param filepath:
        :param filetype: Format of consensus file (accepted by Biopython) default: FASTA
        :return: records
        """
        logging.debug("Loading consensus %s %s" % (filetype, filepath))
        if not os.path.exists(filepath):
            logging.error("Consensus filepath %s does not exist!" % filepath)
        records = list(SeqIO.parse(filepath, filetype))
        logging.debug("The consensus file contains %d records" % len(records))
        assert (len(records) > 0)
        if len(records) > 0:
            logging.warning("At present, genofunk only considers the first sequence in the consensus file")
        return records

    def apply_loaded_edits(self):
        if len(self.edits.edits) == 0:
            return
        for edit in self.edits.edits.reverse():
            record = self.consensus_sequence[edit.sequence_id]
            edit.apply_edit(record)

    def load_input_files(self, reference_filepath, consensus_filepath, edit_filepath=""):
        """
        Load reference JSON, consensus FASTA and edit file (if it exists) and store them
        :param reference_filepath:
        :param consensus_filepath:
        :param edit_filepath:
        :return:
        """
        self.consensus_sequence = self.load_consensus_sequence(consensus_filepath)
        self.reference_info = self.load_reference_info(reference_filepath)
        self.edits = EditFile(edit_filepath)
        self.apply_loaded_edits()