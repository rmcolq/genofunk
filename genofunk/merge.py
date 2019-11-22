import os
import logging
import glob
import json
from Bio import SeqIO

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

class Merge():
    def __init__(self):
        self.consensus_sequence = None
        self.edits = None
        self.coordinates=None

    def load_coordinates_file(self, coordinates_file, features_list=None):
        if not os.path.exists(coordinates_file):
            logging.error("Paired coordinates file %s does not exist!" % coordinates_file)

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

    def load_edits_in_range(self, edit_file, features_list=None):
        assert self.coordinates

        if not features_list:
            self.edits.append(edit_file)
            return

        new_edits = EditFile(edit_file)
        for edit in new_edits.edits:
            for feature in self.coordinates:
                if edit.sequence_id in self.coordinates[feature] and \
                        self.coordinates[feature][edit.sequence_id]["start"] <= edit.sequence_position <= \
                        self.coordinates[feature][edit.sequence_id]["end"]:
                    self.edits.add_edit(edit)

    def load_from_directory(self, directory, features_list=None, filetype="fasta"):
        """
        Looks for pairs of *.fasta, *.fasta.edit files in a directory, and loads them into a single big list of edits
        and a single big list of consensus sequences
        :param directory:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
        self.coordinates = {}
        self.consensus_sequence = []
        self.edits = EditFile()

        edit_files = glob.glob("%s/*.edits" %directory)
        if len(edit_files) == 0:
            logging.error("No edit files found in directory %s" %directory)
            assert(len(edit_files) > 0)
        for edit_file in edit_files:
            consensus_file = edit_file.replace(".edits","")
            if not os.path.exists(consensus_file):
                logging.error("Paired consensus file %s does not exist!" % consensus_file)
                assert(os.path.exists(consensus_file))
            coordinates_file = edit_file.replace(".edits", ".coordinates")

            logging.debug("Loading consensus file %s, edit file %s and coordinates file %s" %(consensus_file, edit_file,
                                                                                              coordinates_file))
            self.consensus_sequence.extend(list(SeqIO.parse(consensus_file, filetype)))
            self.load_coordinates_file(coordinates_file, features_list)
            self.load_edits_in_range(edit_file, features_list)
            logging.debug("Now have %d consensus records and %d edits" %(len(self.consensus_sequence),
                                                                         len(self.edits.edits)))

        self.edits.sort()

    def query_edit(self, edit):
        """
        Asks the user if they wish to ignore or accept a proposed edit.
        :param edit:
        :return:
        """
        response = False
        data = input("Do you wish to ignore edit %s? [Y/n]" %edit)
        while not response:
            if data in ('\n', 'Y', 'y'):
                edit.edit_accepted = False
                response = True
            elif data in ('N', 'n'):
                edit.edit_accepted = True
                response = True
            else:
                data = input("Not an appropriate choice. Do you wish to ignore edit? [Y/n]")

    def query_edits(self, edit_list, message=""):
        """
        Asks the user if they wish to ignore or accept an edit which occurs between a position in a reference and
        several consensus sequences. Applies the judgement to all these edits.
        :param edit_list:
        :param message:
        :return:
        """
        response = False
        data = input("%s\nDo you wish to ignore them? [Y/n]" %message)
        while not response:
            if data in ('\n', 'Y', 'y'):
                for edit in edit_list:
                    edit.edit_accepted = False
                response = True
            elif data in ('N', 'n'):
                for edit in edit_list:
                    edit.edit_accepted = True
                response = True
            else:
                data = input("Not an appropriate choice. Do you wish to ignore them? [Y/n]")

    def add_query_to_edits(self, edit_list):
        for edit in edit_list:
            edit.edit_query = True

    def find_common_edits(self, min_occurrence=2, interactive=False):
        """
        Search for edits which occur with respect to the same reference at the same position and are of the same from/to
        form in multiple consensus records. Questions the user about accepting or ignoring each such case.
        :param min_occurrence: minimum number of times to see an edit to question it (Default:2)
        :return:
        """
        self.edits.sort()
        current_identical = []
        for new_edit in self.edits.edits:
            if len(current_identical) == 0 or (new_edit.reference_id == current_identical[-1].reference_id
                                     and new_edit.reference_position == current_identical[-1].reference_position
                                     and new_edit.edit_from == current_identical[-1].edit_from
                                     and new_edit.edit_to == current_identical[-1].edit_to):
                current_identical.append(new_edit)
            else:
                if len(current_identical) >= min_occurrence:
                    message = "Found %d identical edits like %s.\nThese are likely to be real (not " \
                              "caused by sequencing/assembly errors)." % (len(current_identical), current_identical[-1])
                    if interactive:
                        self.query_edits(current_identical, message)
                    else:
                        self.add_query_to_edits(current_identical)
                current_identical = [new_edit]

    def find_similar_edits(self, min_occurrence=2, interactive=False):
        """
        Search for edits which occur with respect to the same reference at the same position but possibly with DIFFERENT
        to/from sequences in multiple consensus records. Questions the user about accepting or ignoring each such case.
        :param min_occurrence: minimum number of times to see an edit to question it (Default:2)
        :return:
        """
        self.edits.sort()
        current_similar = []
        for new_edit in self.edits.edits:
            if len(current_similar) == 0 or (new_edit.reference_id == current_similar[-1].reference_id
                                             and new_edit.reference_position == current_similar[-1].reference_position):
                current_similar.append(new_edit)
            else:
                if len(current_similar) >= min_occurrence:
                    if interactive:
                        print("Found %d similar edits:" %len(current_similar))
                        for edit in current_similar:
                            print(edit)
                        print("These may be real or caused by sequencing/assembly errors")
                        for edit in current_similar:
                            if edit.edit_accepted:
                                    self.query_edit(edit)
                    else:
                        self.add_query_to_edits(current_similar)
                current_similar = [new_edit]

    def run(self, directory, output_file="all.edits", features="", min_occurence=2, interactive=False):
        if features:
            feature_list = features.split(",")
        else:
            feature_list = None
        self.load_from_directory(directory, feature_list)
        self.find_common_edits(min_occurrence=min_occurence, interactive=interactive)
        self.find_similar_edits(min_occurrence=min_occurence, interactive=interactive)
        if self.edits.contains_query():
            self.edits.save("tmp.%s" % output_file, filter_by_applied=False)
        else:
            self.edits.save(output_file, filter_by_applied=False)