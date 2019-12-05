import os
import sys
import logging
import glob
import json
from Bio import SeqIO
from Bio.Seq import Seq

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

class Apply():
    def __init__(self):
        self.consensus_sequence = None
        self.edits = None
        self.coordinates = None

    def load_coordinates_file(self, coordinates_file, features_list=None):
        """
        Loads coordinates file generated by `annotate` to self.coordinates. If a features_list is provided, only loads
        coordinates that lie within these features
        :param coordinates_file: file of coordinates output by `annotate`
        :param features_list: list of features named in reference JSON to restrict to
        :return:
        """
        assert self.coordinates is not None

        if not os.path.exists(coordinates_file):
            logging.error("Paired coordinates file %s does not exist!" % coordinates_file)
            assert (os.path.exists(coordinates_file))

        if features_list:
            for key in features_list:
                if key not in self.coordinates:
                    self.coordinates[key] = {}

        with open(coordinates_file) as json_file:
            data = json.load(json_file)
            for key in data:
                if features_list and key not in features_list:
                    continue
                if key not in self.coordinates:
                    self.coordinates[key] = {}
                for record_name in data[key]:
                    if record_name in self.coordinates[key]:
                        logging.error(
                            "Record name %s exists in multiple consensus files (and coordinates files). This will "
                            "break things!" % record_name)
                    assert record_name not in self.coordinates[key]
                self.coordinates[key].update(data[key])

    def load_consensus_file(self, consensus_file, filetype="fasta"):

        assert self.consensus_sequence is not None

        if not os.path.exists(consensus_file):
            logging.error("Paired consensus file %s does not exist!" % consensus_file)
            assert (os.path.exists(consensus_file))

        record_dict = SeqIO.index(consensus_file, filetype)
        assert len(record_dict) > 0
        for key in record_dict:
            if key in self.consensus_sequence:
                logging.error(
                    "Record name %s exists in multiple consensus files. This will "
                    "break things!" % key)
            assert key not in self.consensus_sequence
        self.consensus_sequence.update(record_dict)
        logging.debug("After loading consensus file %s, we have %d consensus records " % (consensus_file,
                                                                                          len(self.consensus_sequence)))

    def load_edits_in_range(self, edit_file, features_list=None):
        """
        Loads list of edits from an editfile and appends them to self.edits. If a features_list is provided, only loads
        edits that lie within the coordinate range for those features
        :param edit_file:
        :param features_list: list of features named in reference JSON to restrict to
        :return:
        """
        assert self.coordinates is not None

        if not features_list:
            self.edits = EditFile(edit_file)
        else:
            self.edits = EditFile()
            new_edits = EditFile(edit_file)
            for edit in new_edits.edits:
                for feature in self.coordinates:
                    if edit.sequence_id in self.coordinates[feature] and \
                            self.coordinates[feature][edit.sequence_id]["start"] <= edit.sequence_position <= \
                            self.coordinates[feature][edit.sequence_id]["end"]:
                        self.edits.add_edit(edit)

        self.edits.sort(reverse=True, seq_position=True)
        logging.debug("Now have %d sorted edits" % len(self.edits.edits))

    def load_input_files(self, directory, edit_filepath, features_list=None, filetype="fasta"):
        """
        Looks for pairs of *.fasta, *.fasta.edit, *.fasta.coordinates files in a directory, and loads them into single
        big lists of consensus_sequences and coordinates
        :param directory:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
        self.coordinates = {}
        self.consensus_sequence = {}

        edit_files = glob.glob("%s/*.edits" %directory)
        edit_files = list(filter(lambda x: not x.endswith(edit_filepath), edit_files))

        if len(edit_files) == 0:
            logging.error("No edit files found in directory %s" %directory)
            assert(len(edit_files) > 0)

        for edit_file in edit_files:
            coordinates_file = edit_file.replace(".edits", ".coordinates")
            self.load_coordinates_file(coordinates_file, features_list=features_list)
            consensus_file = edit_file.replace(".edits","")
            self.load_consensus_file(consensus_file, filetype)

        self.load_edits_in_range(edit_filepath, features_list=features_list)

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

    def save_updated_consensuses(self, filepath=None, feature=None, amino_acid=False):
        if filepath:
            out_handle = open(filepath, 'w')
        else:
            out_handle = sys.stdout
            if feature:
                out_handle.write(feature)

        for seq_name in self.consensus_sequence:
            seq = self.consensus_sequence[seq_name]
            assert type(seq), Seq
            if feature:
                if seq_name not in self.coordinates[feature]:
                    continue
                (start, end) = self.coordinates[feature][seq_name]["start"], self.coordinates[feature][seq_name]["end"]
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

    def run(self, directory, edit_filepath, output_prefix, features=""):
        if features:
            features_list = features.split(",")
        else:
            features_list = None
        self.load_input_files(directory, edit_filepath, features_list=features_list)
        self.apply_loaded_edits()

        if features_list:
            for feature in features_list:
                self.save_updated_consensuses("%s.na.fasta" % output_prefix, feature=feature)
                self.save_updated_consensuses("%s.aa.fasta" % output_prefix, feature=feature, amino_acid=True)
        else:
            self.save_updated_consensuses("%s.na.fasta" % output_prefix)
            self.save_updated_consensuses("%s.aa.fasta" % output_prefix, amino_acid=True)
