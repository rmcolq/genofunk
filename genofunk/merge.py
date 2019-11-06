import os
import logging
import glob
from Bio import SeqIO

from genofunk import editfile

EditFile = editfile.EditFile
Edit = editfile.Edit

class Merge():
    def __init__(self):
        self.consensus_sequence = None
        self.edits = None

    def load_from_directory(self, directory, filetype="fasta"):
        """
        Looks for pairs of *.fasta, *.fasta.edit files in a directory, and loads them into a single big list of edits
        and a single big list of consensus sequences
        :param directory:
        :param filetype: if consensus sequence file not FASTA
        :return:
        """
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
            logging.debug("Loading consensus file %s and edit file %s" %(consensus_file, edit_file))
            self.consensus_sequence.extend(list(SeqIO.parse(consensus_file, filetype)))
            self.edits.append(edit_file)
            logging.debug("Now have %d consensus records and %d edits" %(len(self.consensus_sequence), len(self.edits.edits)))

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

    def find_common_edits(self, min_occurrence=2):
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
                    self.query_edits(current_identical, message)
                current_identical = [new_edit]

    def find_similar_edits(self, min_occurrence=2):
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
                    print("Found %d similar edits:" %len(current_similar))
                    for edit in current_similar:
                        print(edit)
                    print("These may be real or caused by sequencing/assembly errors")
                    for edit in current_similar:
                        if edit.edit_accepted:
                            self.query_edit(edit)
                current_similar = [new_edit]

    def run(self, directory):
        self.load_from_directory(directory)
        self.find_common_edits()
        self.find_similar_edits()
        self.edits.save("tmp.edits", filter_by_applied=False)