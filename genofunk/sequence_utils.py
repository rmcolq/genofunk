import logging

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

def str_coordinates(coordinates):
    if coordinates is None:
        return coordinates
    else:
        return ",".join([str(i) for i in coordinates])

def is_open_reading_frame(amino_acid_sequence, allow_stop_codons=False, allow_missing=True, stop_codons=['*']):
    if len(amino_acid_sequence) < 3:
        return False

    if not allow_stop_codons:
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
