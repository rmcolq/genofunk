import logging
from Bio import SeqIO
from Bio.Seq import Seq


def add_key_to_coordinate_dict(coordinate_dict, key):
    if key in coordinate_dict:
        return
    else:
        coordinate_dict[key] = {}

def get_coordinates_from_json(json_entry, pairs=True):
    coordinates = []
    if 'join' in json_entry:
        for coords in json_entry['join']:
            if pairs:
                coordinates.append([coords['start'], coords['end']])
            else:
                coordinates.append(coords['start'])
                coordinates.append(coords['end'])
    else:
        if pairs:
            coordinates.append([json_entry['start'], json_entry['end']])
        else:
            coordinates.append(json_entry['start'])
            coordinates.append(json_entry['end'])
    return coordinates

def apply_edits_in_range(editfile, record, offset=0, coordinates=None, filter_by_accepted=True):
    """
    Apply the edits to the consensus nucleotide sequences in place (in reverse order to avoid offset errors)
    :return:
    """
    if record is None or editfile is None or len(editfile.edits) == 0:
        return record, 0
    coordinate_difference = 0
    if len(editfile.edits) == 0:
        return record
    for edit in editfile.edits:
        if edit.sequence_id == record.id and \
                (coordinates is None or coordinates[0] <= edit.sequence_position < coordinates[-1]):
            record, applied = edit.apply_edit(record, offset=offset+coordinate_difference,
                                              filter_by_accepted=filter_by_accepted)
            if applied:
                coordinate_difference += len(edit.edit_to) - len(edit.edit_from)
            logging.debug("coordinate_difference has been updated by %i and is now %i" %(len(edit.edit_to) -
                                                                                         len(edit.edit_from),
                                                                                         coordinate_difference))
    return record, coordinate_difference

def make_sequence_length_divide_by_3(sequence):
    if len(sequence) % 3 == 1:
        sequence = sequence + "NN"
    elif len(sequence) % 3 == 2:
        sequence = sequence + "N"
    return sequence

def count_stops(sequence, stop_codon="*"):
    return sequence.count(stop_codon)

def find_run_n( sequence, min_run_length=9):
    i = sequence.find("N")
    j = 0
    while i >= 0:
        for j in range(i + 1, len(sequence)):
            if sequence[j] != "N":
                break
        if j - 1 - i >= min_run_length:
            return i,j
        i = sequence.find("N", j)
    return -1, -1

def rfind_run_n( sequence, min_run_length=9):
    i = sequence.rfind("N")
    j = 0
    while i >= 0:
        for j in reversed(range(0, i)):
            if sequence[j] != "N":
                break
        if i - j - 1 >= min_run_length:
            return j + 1, i+1
        i = sequence.rfind("N", 0, j)
    return -1, -1

def find_n_runs(sequence, min_run_length=3):
    runs = []
    run_start = None
    run_end = None
    for i,c in enumerate(sequence):
        if c in ["X", "N"] and run_start is None:
            run_start = i
            run_end = i
        elif c in ["X", "N"]:
            run_end += 1
        elif run_start is not None:
            if run_end - run_start >= min_run_length:
                runs.append([run_start, run_end])
            run_start = None
            run_end = None
    return runs

def shift_nucleotide_sequence_into_frame(sequence, coordinates=None):
    num_stop_codons = []

    shift_0 = make_sequence_length_divide_by_3(sequence)
    translated_0 = shift_0.translate()
    num_stop_codons.append(count_stops(translated_0))

    shift_1 = "N" + sequence
    shift_1 = make_sequence_length_divide_by_3(shift_1)
    translated_1 = shift_1.translate()
    num_stop_codons.append(count_stops(translated_1))

    shift_2 = "NN" + sequence
    shift_2 = make_sequence_length_divide_by_3(shift_2)
    translated_2 = shift_2.translate()
    num_stop_codons.append(count_stops(translated_2))

    seqs = [shift_0, shift_1, shift_2]
    i = num_stop_codons.index(min(num_stop_codons))
    if coordinates is not None:
        coordinates = [coordinates[0] - i, coordinates[1] - i]
    return seqs[i], coordinates

def codon_aware_update(old_coordinate, new_coordinate):
    while (new_coordinate - old_coordinate) % 3 != 0:
        new_coordinate += 1
    return new_coordinate

def get_sequence(seq, coordinates=None, shift_into_frame=False, offset=None, amino_acid=True, stop_symbol='*'):
    """
    Takes a Biopython Seq object and subsets the sequence within a coordinate range, subsequently offsetting and
    translating into amino acid sequence as required
    :param seq:
    :param coordinates: 0-based (start,end) in nucleotide sequence (end NOT included)
    :param offset:
    :param amino_acid:
    :return: string
    """
    assert(type(seq), Seq)

    if coordinates:
        start, end = coordinates[0], coordinates[-1]
        while start < 0:
            seq = "N" + seq
            start += 1
        end = min(end, len(seq))

        return_sequence = ""
        updated_coordinates = [start]
        if len(coordinates) > 2:
            updated_coordinates.extend(coordinates[1:-1])
        updated_coordinates.append(end)
        #logging.debug("updated coordinates %s" %str_coordinates(updated_coordinates))
        i = 0
        while i < len(updated_coordinates):
            (start, end) = updated_coordinates[i:i + 2]
            return_sequence += seq[start:end]
            i = i + 2
        seq = return_sequence
    if offset:
        seq = seq[offset:]
    if shift_into_frame:
        seq, coordinates = shift_nucleotide_sequence_into_frame(seq, coordinates)
    else:
        seq = make_sequence_length_divide_by_3(seq)
    if "-" in str(seq):
        seq = Seq(str(seq).replace("-","N"))
    if amino_acid:
        seq = seq.translate(stop_symbol=stop_symbol)
    return str(seq), coordinates

def str_coordinates(coordinates):
    if coordinates is None:
        return coordinates
    else:
        return ",".join([str(i) for i in coordinates])

def get_stop_codon_positions(query_sequence, stop_codons, min_frame_shift_position):
    positions = []
    for stop in stop_codons:
        position = query_sequence.find(stop, min_frame_shift_position)
        logging.debug("searching for %s in %s after position %i" % (stop, query_sequence, min_frame_shift_position))
        if position > 0:
            positions.append(position)
            logging.debug("Found stop position %i" % position)
    return positions

def get_num_stop_codons(query_sequence, stop_codons, min_frame_shift_position):
    count = 0
    for stop in stop_codons:
        count += query_sequence.count(stop, min_frame_shift_position)
    return count

def has_fewer_stop_codons(old_sequence, new_sequence, stop_codons):
    old_count = get_num_stop_codons(old_sequence[:-2], stop_codons, min_frame_shift_position=0)
    new_count = get_num_stop_codons(new_sequence[:-2], stop_codons, min_frame_shift_position=0)
    if new_count <= old_count:
        return True
    return False

def is_open_reading_frame(amino_acid_sequence, allow_stop_codons_in_middle=False, allow_missing=True,
                          stop_codons=['*']):
    if len(amino_acid_sequence) < 3:
        return False

    if not allow_stop_codons_in_middle:
        for letter in amino_acid_sequence[1:-1]:
            if letter in stop_codons:
                return False

    start = ["M"]
    stop = stop_codons[:]
    if allow_missing:
        start.append("X")
        stop.append("X")

    if amino_acid_sequence[0] not in start:
        return False

    if amino_acid_sequence[-1] not in stop:
        return False

    return True
