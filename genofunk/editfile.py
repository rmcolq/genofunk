import pandas as pd
import logging
import functools

@functools.total_ordering
class Edit():
    def __init__(self, sequence_id, sequence_position, edit_from, edit_to, reference_id, reference_position, edit_accepted=True, edit_applied=False, edit_query=False):
        self.sequence_id = sequence_id
        self.sequence_position = sequence_position
        self.edit_from = edit_from
        self.edit_to = edit_to
        self.reference_id = reference_id
        self.reference_position = reference_position
        self.edit_accepted = edit_accepted
        self.edit_applied = edit_applied
        self.edit_query = edit_query
        
    def __repr__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def __eq__(self, other):
        if self.sequence_id == other.sequence_id and self.sequence_position == other.sequence_position \
                and self.edit_from == other.edit_from and self.edit_to == other.edit_to \
                and self.reference_id == other.reference_id and self.reference_position == other.reference_position:
            return True
        return False

    def __ne__(self, other):
        if self.sequence_id == other.sequence_id and self.sequence_position == other.sequence_position \
                and self.edit_from == other.edit_from and self.edit_to == other.edit_to \
                and self.reference_id == other.reference_id and self.reference_position == other.reference_position:
            return False
        return True

    def __lt__(self, other):
        if self.reference_id < other.reference_id:
            return True
        elif self.reference_id > other.reference_id:
            return False
        elif self.reference_position < other.reference_position:
            return True
        elif self.reference_position > other.reference_position:
            return False
        elif self.edit_from < other.edit_from:
            return True
        elif self.edit_from > other.edit_from:
            return False
        elif self.edit_to < other.edit_to:
            return True
        elif self.edit_to > other.edit_to:
            return False
        elif self.sequence_id < other.sequence_id:
            return True
        elif self.sequence_id > other.sequence_id:
            return False
        elif self.sequence_position < other.sequence_position:
            return True
        elif self.sequence_position > other.sequence_position:
            return False
        return False

    def apply_edit(self, record, offset=0, filter_by_accepted=False):
        """
        Applies the edit to the consensus nucleotide sequence in place.
        :param record: a Seq object
        :param offset: account for previously applied edits earlier in the sequence which
        :return:
        """
        print("apply edit", self)
        print(self.edit_applied)
        if self.edit_applied:
            return

        print(filter_by_accepted, self.edit_accepted)
        if filter_by_accepted and not self.edit_accepted:
            return

        logging.debug(record.id)
        sequence = record.seq
        if self.edit_from == "N":
            self.edit_from = sequence[self.sequence_position + offset]
        elif self.edit_from == "NN":
            self.edit_from = str(sequence[self.sequence_position + offset:self.sequence_position + offset + 2])

        assert(sequence[self.sequence_position + offset:self.sequence_position + offset + len(self.edit_from)] == self.edit_from)
        updated_sequence = sequence[:self.sequence_position + offset]
        updated_sequence += self.edit_to 
        updated_sequence += sequence[self.sequence_position + offset + len(self.edit_from):]
    
        record.seq = updated_sequence
        self.edit_applied = True
        logging.debug("Applying edit changes")
        logging.debug("Old: %s" %sequence)
        logging.debug("New: %s" %updated_sequence)
        
    def remove_edit(self, record, offset=0):
        """
        Removes the edit to the consensus nucleotide sequence in place.
        :param record: a Seq object
        :param offset: account for previously applied edits earlier in the sequence which
        :return:
        """
        if not self.edit_applied:
            return

        sequence = record.seq
        assert(sequence[self.sequence_position + offset:self.sequence_position + offset + len(self.edit_to)]
               == self.edit_to)
        updated_sequence = sequence[:self.sequence_position + offset]
        updated_sequence += self.edit_from 
        updated_sequence += sequence[self.sequence_position + offset + len(self.edit_to):]
    
        record.seq = updated_sequence
        self.edit_applied = False
        
class EditFile():
    def __init__(self, filepath=None):
        self.edits = []
        if filepath:
            self.append(filepath)
    
    def __repr__(self):
        return  str(self.__class__) + '\n'+ '\n'.join(('{} = {}'.format(item, self.__dict__[item]) for item in
                                                       self.__dict__))
    def __eq__(self, other):
        if len(self.edits) != len(other.edits):
            return False
        self.sort()
        other.sort()
        for i in range(len(self.edits)):
            if self.edits[i] != other.edits[i]:
                return False
        return True

    def add_edit(self, edit):
        """
        Add edit to list for file
        :param edit:
        :return:
        """
        if edit not in self.edits:
            self.edits.append(edit)

    def append(self, filepath):
        """
        Append edits found in filepath to list of edits
        :param filepath:
        :return:
        """
        data = pd.read_csv(filepath)
        for i,row in data.iterrows():
            edit_from, edit_to = str(row["from"]), str(row["to"])
            if edit_from in [None, "nan"]:
                edit_from = ""
            if edit_to in [None, "nan"]:
                edit_to = ""
            edit_accepted = True
            if 'edit_accepted' in data.columns:
                edit_accepted = row['edit_accepted']
            edit_query = False
            if 'query' in data.columns:
                edit_query = row['query']

            e = Edit(row["read_id"], row["read_pos"], edit_from, edit_to, row["ref_id"], row["ref_pos"],
                     edit_accepted=edit_accepted, edit_query=edit_query)
            self.edits.append(e)

    def sort(self, reverse=False):
        """
        Sorts edits by reference_id, reference_position, edit_from, edit_to, sequence_id, sequence_position
        :return:
        """
        self.edits.sort(reverse=reverse)

    def contains_query(self):
        for edit in self.edits:
            if edit.edit_query:
                return True
        return False

    def save(self, filepath, filter_by_applied=True, filter_by_accepted=False):
        """
        Save EditFile as a CSV
        :param filepath:
        :param filter_by_applied: filter out edits which have not been applied? Default: True
        :return:
        """
        columns = ["read_id", "read_pos", "from", "to", "ref_id", "ref_pos", "edit_accepted"]
        if self.contains_query():
            columns.append("query")
        with open(filepath, "w") as f:
            header = ','.join(columns)
            f.write("%s\n" %header)
            for e in self.edits:
                if (not filter_by_applied or e.edit_applied) and (not filter_by_accepted or e.edit_accepted):
                    attributes = ','.join(
                        [str(e.sequence_id), str(e.sequence_position), e.edit_from, e.edit_to, e.reference_id,
                         str(e.reference_position), str(e.edit_accepted)])
                    if self.contains_query() and e.edit_query:
                        attributes += ",?"
                    elif self.contains_query():
                        attributes += ","
                    f.write("%s\n" %attributes)