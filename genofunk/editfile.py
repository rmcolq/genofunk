import pandas as pd
import logging


class Edit():
    def __init__(self, sequence_id, sequence_position, edit_from, edit_to, reference_id, reference_position):
        self.sequence_id = sequence_id
        self.sequence_position = sequence_position
        self.edit_from = edit_from
        self.edit_to = edit_to
        self.reference_id = reference_id
        self.reference_position = reference_position
        self.edit_applied = False
        
    def __repr__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def apply_edit(self, record, offset=0):
        sequence = record.seq
        if self.edit_from == "N":
            self.edit_from = sequence[self.sequence_position + offset]
        elif self.edit_from == "NN":
            self.edit_from = str(sequence[self.sequence_position + offset:self.sequence_position + offset + 2])

        assert(sequence[self.sequence_position + offset:self.sequence_position + offset + len(self.edit_from)]
               == self.edit_from)
        updated_sequence = sequence[:self.sequence_position + offset]
        updated_sequence += self.edit_to 
        updated_sequence += sequence[self.sequence_position + offset + len(self.edit_from):]
    
        record.seq = updated_sequence
        self.edit_applied = True
        logging.debug("Applying edit changes")
        logging.debug("Old: %s" %sequence)
        logging.debug("New: %s" %updated_sequence)
        
    def remove_edit(self, record, offset=0):
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
            data = pd.read_csv(filepath)
            for i,row in data.iterrows():
                e = Edit(row["read_id"], row["read_pos"], row["from"], row["to"], row["ref_id"], row["ref_pos"])
                self.edits.append(e)
    
    def __repr__(self):
        return  str(self.__class__) + '\n'+ '\n'.join(('{} = {}'.format(item, self.__dict__[item]) for item in
                                                       self.__dict__))
            
    def add_edit(self, edit):
        self.edits.append(edit)

    def save(self, filepath):
        with open(filepath, "w") as f:
            header = ','.join(["read_id", "read_pos", "from", "to", "ref_id", "ref_pos"])
            f.write("%s\n" %header)
            for e in self.edits:
                if e.edit_applied:
                    attributes = ','.join(
                        [str(e.sequence_id), str(e.sequence_position), e.edit_from, e.edit_to, e.reference_id,
                         str(e.reference_position)])
                    f.write("%s\n" %attributes)