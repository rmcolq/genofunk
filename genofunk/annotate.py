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


class Annotate:
    def __init__(self, closest_accession=None):
        self.consensus_sequence = None
        self.reference_info = None
        self.edits = None
        self.coordinates = {}
        self.closest_accession = closest_accession
        
    def load_reference_info(self, filepath):
        """
        Load the reference JSON, checking it contains the appropriate fields
        :param filepath:
        :return: dict (from JSON)
        """
        logging.debug("Loading reference JSON %s" % filepath)
        if not filepath or not os.path.exists(filepath):
            logging.error("Reference filepath %s does not exist!" %filepath)
            assert filepath
            assert (os.path.exists(filepath))
        with open(filepath) as json_file:
            data = json.load(json_file)
        logging.debug("Checking that JSON has correct format and contains the appropriate fields"
                      "for accession %s" %self.closest_accession)
        assert ('references' in data.keys())
        assert ('features' in data.keys())
        assert (len(data['features']) > 0)
        if self.closest_accession is None and len(data['references'].keys()) == 1:
            self.closest_accession = list(data['references'].keys())[0]
        assert (self.closest_accession in data['references'].keys())
        assert ('sequence' in data['references'][self.closest_accession].keys())
        assert ('locations' in data['references'][self.closest_accession].keys())
        assert (len(data['references'][self.closest_accession]['locations']) > 0)
        return data
        
    def load_consensus_sequence(self, filepath, filetype="fasta"):
        """
        Load consensus FASTA (or other file type) and check there is at least one record
        :param filepath:
        :param filetype: Format of consensus file (accepted by Biopython) default: FASTA
        :return: records
        """
        logging.debug("Loading consensus %s %s" %(filetype,filepath))
        if not filepath or not os.path.exists(filepath):
            logging.error("Consensus filepath %s does not exist!" %filepath)
            assert (filepath != None)
            assert(os.path.exists(filepath))
        records = list(SeqIO.parse(filepath, filetype))
        logging.debug("The consensus file contains %d records" %len(records))
        assert(len(records) > 0)
        #if len(records) > 0:
        #    logging.warning("At present, genofunk only considers the first sequence in the consensus file")
        return records

    def apply_loaded_edits(self):
        """
        Apply the edits to the consensus nucleotide sequences in place (in reverse order to avoid offset errors)
        :return:
        """
        if len(self.edits.edits) == 0:
            return
        for edit in self.edits.edits.reverse():
            record = self.consensus_sequence[edit.sequence_id]
            edit.apply_edit(record)
    
    def load_input_files(self, reference_filepath, consensus_filepath, edit_filepath = ""):
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

    def make_sequence_length_divide_by_3(self, sequence):
        if len(sequence) % 3 == 1:
            sequence = sequence + "NN"
        elif len(sequence) % 3 == 2:
            sequence = sequence + "N"
        return sequence

    def count_stops(self, sequence, stop_codon="*"):
        return sequence.count(stop_codon)

    def shift_nucleotide_sequence_into_frame(self, sequence, coordinates=None):
        num_stop_codons = []

        shift_0 = self.make_sequence_length_divide_by_3(sequence)
        translated_0 = shift_0.translate()
        num_stop_codons.append(self.count_stops(translated_0))

        shift_1 = "N" + sequence
        shift_1 = self.make_sequence_length_divide_by_3(shift_1)
        translated_1 = shift_1.translate()
        num_stop_codons.append(self.count_stops(translated_1))

        shift_2 = "NN" + sequence
        shift_2 = self.make_sequence_length_divide_by_3(shift_2)
        translated_2 = shift_2.translate()
        num_stop_codons.append(self.count_stops(translated_2))

        seqs = [shift_0, shift_1, shift_2]
        i = num_stop_codons.index(min(num_stop_codons))
        if coordinates is not None:
            coordinates = [coordinates[0] - i, coordinates[1] - i]
        return seqs[i], coordinates


    def get_sequence(self, seq, coordinates=None, shift_into_frame=False, offset=None, amino_acid=True):
        """
        Takes a Biopython Seq object and subsets the sequence within a coordinate range, subsequently offsetting and
        translating into amino acid sequence as required
        :param seq:
        :param coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param offset:
        :param amino_acid:
        :return: string
        """
        assert(type(seq), Seq)
        if coordinates:
            (start, end) = coordinates
            while start < 0:
                seq = "N" + seq
                start += 1
            seq = seq[start:end]
        if offset:
            seq = seq[offset:]
        if shift_into_frame:
            seq, coordinates = self.shift_nucleotide_sequence_into_frame(seq, coordinates)
        else:
            seq = self.make_sequence_length_divide_by_3(seq)
        if amino_acid:
            seq = seq.translate()
        return str(seq), coordinates
    
    def get_reference_sequence(self, coordinates=None, shift_into_frame=False, offset=None, amino_acid=True, accession=None):
        """
        Get the string associated to the given/closest accession (within a coordinate range and translated as required)
        :param coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param offset:
        :param amino_acid:
        :param accession:
        :return: string
        """
        if not accession:
            accession = self.closest_accession
            assert(accession)
        seq = Seq(self.reference_info['references'][accession]['sequence'])
        return self.get_sequence(seq, coordinates, shift_into_frame, offset, amino_acid)
    
    def get_query_sequence(self, record_id=0, coordinates=None, shift_into_frame=False, offset=None, amino_acid=True):
        """
        Get the string associated to the given record (within a coordinate range and translated as required)
        :param record_id:
        :param coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param offset:
        :param amino_acid:
        :return:
        """
        seq = self.consensus_sequence[record_id].seq
        return self.get_sequence(seq, coordinates, shift_into_frame, offset, amino_acid)
    
    def decode_cigar(b):
        __BAM_CIGAR_STR = 'MIDNSHP=X8'
        def _decode(x):
            l = str(x>>4)
            try:
                c = __BAM_CIGAR_STR[x&0xf]
            except:
                c = 'M'
            return l+c
        return ''.join([_decode(i) for i in b])
            
    def pairwise_ssw_align(self, ref_seq, query_seq, gap_open=10, gap_extension=1, matrix=parasail.blosum62):
        """
        Run parasail-python SSW function
        :param ref_seq:
        :param query_seq:
        :param gap_open: Default 10
        :param gap_extension: Default 1
        :param matrix: Default BLOSUM62
        :return: parasail result (includes score1, cigar, ref_begin1, ref_end1, read_begin1, read_end1 attributes)
        """
        result = parasail.ssw(query_seq, ref_seq, gap_open, gap_extension, matrix)
        return result 
    
    def pairwise_sw_trace_align(self, ref_seq, query_seq, gap_open=10, gap_extension=1, matrix=parasail.blosum62):
        """

        :param ref_seq:
        :param query_seq:
        :param gap_open: Default 10
        :param gap_extension: Default 1
        :param matrix: Default BLOSUM62
        :return: parasail result (includes a cigar which can be decoded)
        """
        result = parasail.nw_trace(query_seq, ref_seq, gap_open, gap_extension, matrix)
        logging.debug("Parasail result %s" %result.cigar.decode)
        return result 
    
    def parse_cigar(self,result):
        """
        Extract the cigar from the parasail sw_trace_align result and turn it into cigar pairs
        :param result:
        :return: list of pairs of (match/mismatch type, length)
        """
        #M   alignment match (can be a sequence match or mismatch)
        #I   insertion to the reference
        #D   deletion from the reference
        #N   skipped region from the reference
        #S   soft clipping (clipped sequences present in SEQ)
        #H   hard clipping (clipped sequences NOT present in SEQ)
        #P   padding (silent deletion from padded reference)
        #=   sequence match
        #X   sequence mismatch
        pairs = []
        cigar = result.cigar.decode.decode('UTF-8')
        last_position = 0
        current_position = 0
        for c in cigar:
            if not c.isdigit():
                pairs.append((c,int(cigar[last_position:current_position])))
                last_position = current_position+1
            current_position += 1
        return pairs

    def cigar_length(self, pairs, max_mismatch, n_runs=[], min_match=3):
        """
        Find the length of aligned sequence until the first insertion/deletion/padding
        :param pairs:
        :return: number
        """
        #pairs = self.parse_cigar(result)
        total = 0
        subtotal = 0
        found_n_run = False
        have_min_match_after_n_run = False
        for c,i in pairs:
            if c in ["="]:
                total += i + subtotal
                subtotal = 0
                if found_n_run and i >= min_match:
                    have_min_match_after_n_run = True
            elif c in ["M"]:
                subtotal += i
            elif c in ["X"] and i <= max_mismatch and not found_n_run:
                subtotal += i
            elif c in ["X"] and i <= max_mismatch and found_n_run and have_min_match_after_n_run:
                subtotal += i
            else:
                found_n_run = False
                for run in n_runs:
                    if run[0] <= total + subtotal <= run[1]:
                        subtotal += i
                        found_n_run = True
                        have_min_match_after_n_run = False
                        break
                if not found_n_run:
                    break
        return total

    def cigar_has_no_indels(self,pairs):
        """
        Are there non-match/mismatch symbols in cigar?
        :param pairs:
        :return: bool
        """
        total = 0
        for c,i in pairs:
            if c in ["I", "D", "N", "S", "H", "P"]:
                return False
        return True

    def cigar_score(self,pairs, match_score=1, mismatch_score=-1, gap_score=-1):
        """
        Find a custom score for cigar
        :param pairs:
        :return: number
        """
        total = 0
        for c,i in pairs:
            if c in ["="]:
                total += i*match_score
            elif c in ["X"]:
                total += i * mismatch_score
            elif c in ["I", "D", "N", "S", "H", "P"]:
                total += i * gap_score
        return total

    def case(self, c):
        #helper function
        if c in ["I", "D", "N", "S", "H", "P"]:
            return 1
        elif c in ["M", "X"]:
            return 2
        elif c == "=":
            return 3
        else:
            return 0

    def is_extended_cigar_prefix(self, old_cigar_pairs, new_cigar_pairs):
        """
        Checks if new_cigar_pairs is just an extension of old_cigar_pairs e.g. ("=",5) -> ("=",5),("X",1),("=",2)
        or ("=",5) -> ("=",6).
        :param old_cigar_pairs:
        :param new_cigar_pairs:
        :return: bool
        """
        if len(old_cigar_pairs) > len(new_cigar_pairs):
            return False

        for i, (old_c, old_c_length) in enumerate(old_cigar_pairs):
            if new_cigar_pairs[i] == old_cigar_pairs[i]:
                continue
            elif i == len(old_cigar_pairs):
                new_c, new_c_length = new_cigar_pairs[i]
                if new_c == old_c and old_c_length < new_c_length:
                    return True
            else:
                return False
        return len(old_cigar_pairs) < len(new_cigar_pairs)

    def is_improved_cigar_prefix(self, old_cigar_pairs, new_cigar_pairs):
        """
        Checks whether the new cigar is thought to be an improvement over the old one. This is determined by checking if
        the new cigar extends the old one (in which case it is an improvement), whether one has a frame shift and the
        other doesn't (frame shifts are bad) and then carefully pairwise checks along cigar prefix to first difference
        and uses that to determine which is better.
        :param old_cigar_pairs:
        :param new_cigar_pairs:
        :return: bool
        """
        logging.debug("Old cigar pairs %s" %old_cigar_pairs)
        logging.debug("New cigar pairs %s" %new_cigar_pairs)

        # First on whether is extension
        if self.is_extended_cigar_prefix(old_cigar_pairs, new_cigar_pairs):
            return True
        elif self.is_extended_cigar_prefix(new_cigar_pairs, old_cigar_pairs):
            return False

        # Second on existance of frameshift
        if self.cigar_has_no_indels(old_cigar_pairs) and not self.cigar_has_no_indels(new_cigar_pairs):
            return False
        elif self.cigar_has_no_indels(new_cigar_pairs) and not self.cigar_has_no_indels(old_cigar_pairs):
            return True

        # Then on prefix comparison
        for i, (old_c, old_c_length) in enumerate(old_cigar_pairs):
            if i >= len(new_cigar_pairs):
                return False
            if new_cigar_pairs[i] == old_cigar_pairs[i]:
                continue

            new_c, new_c_length = new_cigar_pairs[i]

            if new_c == old_c and new_c_length < old_c_length:
                if i + 1 < len(new_cigar_pairs):
                    new_c, new_c_length = new_cigar_pairs[i + 1]
                else:
                    new_c = None
            elif new_c == old_c and new_c_length > old_c_length:
                if i + 1 < len(old_cigar_pairs):
                    old_c, old_c_length = old_cigar_pairs[i + 1]
                else:
                    old_c = None
            logging.debug("Comparing %s and %s and found improvement to be %s" %(new_c, old_c,
                                                                                 self.case(new_c) > self.case(old_c)))
            return self.case(new_c) > self.case(old_c)

        if len(new_cigar_pairs) > len(old_cigar_pairs):
            return True
        return False

    def is_longer_cigar_prefix(self, old_cigar_pairs, new_cigar_pairs, max_mismatch):
        """
        Checks if new cigar prefix is longer than old cigar prefix. Only counts mismatch bases if there is a match after
        and breaks at first indel.
        :param old_cigar_pairs:
        :param new_cigar_pairs:
        :return: bool
        """
        old_cigar_length = self.cigar_length(old_cigar_pairs, max_mismatch)
        new_cigar_length = self.cigar_length(new_cigar_pairs, max_mismatch)
        logging.debug("Comparing %s and %s and found improvement to be %s" % (old_cigar_length, new_cigar_length,
                                                                              old_cigar_length < new_cigar_length))
        return old_cigar_length <= new_cigar_length



        #def choose_closest_reference(self):
    #    closest_accession = None
    #    return closest_accession

    def find_run_n(self, sequence, min_run_length=3):
        i = sequence.find("N")
        j = 0
        while i >= 0:
            for j in range(i + 1, len(sequence)):
                if sequence[j] != "N":
                    break
            if j - 1 - i >= min_run_length:
                return i,j-1
            i = sequence.find("N", j)
        return -1, -1

    def rfind_run_n(self, sequence, min_run_length=3):
        i = sequence.rfind("N")
        j = 0
        while i >= 0:
            for j in reversed(range(0, i)):
                if sequence[j] != "N":
                    break
            if i - j - 1 >= min_run_length:
                return j + 1, i
            i = sequence.rfind("N", 0, j)
        return -1, -1

    def codon_aware_update(self, old_coordinate, new_coordinate):
        while (new_coordinate - old_coordinate) % 3 != 0:
            new_coordinate += 1
        return new_coordinate

    def update_coordinates_if_n_runs(self, query_start, query_end, query_sequence, ref_match_start, ref_match_end,
                                     feature_length):
        first_run_n = self.find_run_n(query_sequence)
        last_run_n = self.rfind_run_n(query_sequence)
        logging.debug("Found first run Ns %i,%i and last run Ns %i,%i" % (first_run_n[0], first_run_n[1], last_run_n[0],
                                                                          last_run_n[1]))

        ref_match_length = ref_match_end - ref_match_start
        tolerance = 32
        if first_run_n[0] == last_run_n[0] \
                and first_run_n[0] < query_end - first_run_n[1] \
                and query_start + first_run_n[0] <= query_end - feature_length + tolerance:
            query_start = self.codon_aware_update(query_start, query_end - feature_length)
            logging.debug("Only one run Ns at start")
        elif first_run_n[0] == last_run_n[0] \
                and first_run_n[0] > query_end - first_run_n[1] \
                and query_start + first_run_n[1] >= query_start + feature_length - tolerance:
            query_end = self.codon_aware_update(query_start, query_start + feature_length)
            logging.debug("Only one run Ns at end")
        elif first_run_n[0] != last_run_n[0] and last_run_n[1] - first_run_n[0] > feature_length:
            query_start = self.codon_aware_update(query_start, query_start + first_run_n[0])
            query_end = self.codon_aware_update(query_start, query_start + last_run_n[1])
            logging.debug("Two runs Ns at start and end")
        else:
            logging.debug("Did not update coordinates")
        return query_start, query_end


    def identify_feature_coordinates(self, feature_coordinates, record_id=0, min_score = 300, max_query_start_offset = 200, max_nucleotide_length_difference = 10):
        """
        Find the region in query/consensus sequence which aligns to sequence in reference nucleotide sequence with
        coordinates feature_coordinates
        :param feature_coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param record_id: Default 0
        :return:
        """
        logging.debug("Had ref feature coordinates %s" % self.str_coordinates(feature_coordinates))
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, shift_into_frame=True,
                                                                        amino_acid=False)
        logging.debug("Updated ref feature coordinates %s" % self.str_coordinates(feature_coordinates))
        offset = max(feature_coordinates[0] - max_query_start_offset, 0)
        query_end = self.codon_aware_update(offset, feature_coordinates[1] + max_query_start_offset)
        query_sequence, coordinates = self.get_query_sequence(record_id,
                                                              coordinates=(offset, query_end),
                                                              amino_acid=False)
        logging.debug("Start with query feature coordinates %s" % self.str_coordinates(coordinates))
        result = self.pairwise_ssw_align(ref_sequence, query_sequence)
        if result and result.score1 > min_score:
            logging.debug("Found score %s and cigar %s" % (result.score1, result.cigar))
            logging.debug("Found query start and end %s, %s" % (result.read_begin1, result.read_end1))
            logging.debug("Found ref start and end %s, %s" % (result.ref_begin1, result.ref_end1))

            query_start = offset + result.read_begin1 - result.ref_begin1
            logging.debug("Update query start to %i" % query_start)
            feature_length = feature_coordinates[1] - feature_coordinates[0]
            query_end = self.codon_aware_update(query_start, offset + result.read_end1 + feature_length - result.ref_end1)
            logging.debug("Update query end to %i" % query_end)
            found_length = query_end - query_start
            query_sequence, coordinates = self.get_query_sequence(record_id,
                                                                  coordinates=(query_start, query_end),
                                                                  amino_acid=False)
            logging.debug("Query sequence %s" %query_sequence)
            logging.debug("Found length %i and feature length %i" %(found_length, feature_length))

            if found_length - feature_length > 0:
                query_start, query_end = self.update_coordinates_if_n_runs(query_start, query_end, query_sequence,
                                                                           result.ref_begin1, result.ref_end1,
                                                                           feature_length)
                logging.debug("Update query start, end to %i, %i" % (query_start, query_end))
                found_length = query_end - query_start
                query_sequence, coordinates = self.get_query_sequence(record_id, coordinates=(query_start, query_end),
                                                                      amino_acid=False)
                query_start, query_end = coordinates
                logging.debug("Query sequence %s" % query_sequence)
                logging.debug("Found length %i and feature length %i" % (found_length, feature_length))

            #if found_length - feature_length > 0 and query_sequence[found_length - feature_length] == "N":
            #    query_start = query_end - feature_length
            #    logging.debug("Update query start to %i to prune off additional bases caused by Ns" % query_start)
            #if (result.read_end1 - result.read_begin1) - (result.ref_end1 - result.ref_begin1) \
            #        > max_nucleotide_length_difference:
            #    query_end = query_start + feature_length - result.ref_begin1
            #    logging.debug("Update query end to %i based on overhang" % query_end)
            return query_start, query_end
        else:
            return None, None

    def find_n_runs(self, sequence, min_run_length=3):
        runs = []
        run_start = None
        run_end = None
        for i,c in enumerate(sequence):
            if c == "X" and run_start is None:
                run_start = i
                run_end = i
            elif c == "X":
                run_end += 1
            elif run_start is not None:
                if run_end - run_start >= min_run_length:
                    runs.append([run_start, run_end])
                run_start = None
                run_end = None
        return runs



    def get_position_for_frame_shift(self, found_coordinates, record_id, cigar_pairs, stop_codons, max_mismatch):
        positions = []
        query_sequence, coordinates = self.get_query_sequence(record_id, coordinates=found_coordinates)
        logging.debug("Have query sequence %s and coordinates %s" % (query_sequence, coordinates))
        for stop in stop_codons:
            position = query_sequence.find(stop)
            if position > 0:
                positions.append(position)
                logging.debug("Found stop codon position %d" % position)
        min_run_length = max_mismatch
        n_runs = self.find_n_runs(query_sequence,min_run_length)
        cigar_length = self.cigar_length(cigar_pairs, max_mismatch, n_runs)
        logging.debug("Found cigar length position %d" % cigar_length)
        positions.append(cigar_length)
        #logging.debug("Found cigar position %d" % positions[-1])
        return min(positions)


    def frame_shift(self, feature_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs, shift_from,
                    shift_to, shift_position, coordinate_difference=0):
        """
        Create a potential edit which applies a frame shift at a given position in the query sequence
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record_id:
        :param ref_sequence: amino acid sequence of reference in features
        :param cigar_pairs: best amino acid alignment of query before frame shift
        :param shift_from:
        :param shift_to:
        :param coordinate_difference: amount to offset when applying frame shift (given that we want raw coordinates in
        query consensus sequence and it may have been updated with previously occurring frame shifts
        :return: updated coordinate_difference, updated cigar_pairs, whether updated, edit (not applied)
        """
        logging.debug("Frame_shift from '%s' to '%s'" %(shift_from, shift_to))

        record_name = self.consensus_sequence[record_id].id
        e = Edit(record_name, found_coordinates[0] + 3 * (shift_position), shift_from, shift_to, self.closest_accession,
                 feature_coordinates[0] + 3 * (shift_position))
        e.apply_edit(self.consensus_sequence[record_id], coordinate_difference)
        updated_coordinate_difference = coordinate_difference + len(shift_to) - len(shift_from)
        updated_found_coordinates = [found_coordinates[0], found_coordinates[1] + updated_coordinate_difference]
        query_sequence, coordinates = self.get_query_sequence(record_id, coordinates=updated_found_coordinates)
        logging.debug("Found query sequence %s and coordinates %s" %(query_sequence, coordinates))
        logging.debug("Have ref sequence %s" % ref_sequence)
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        new_cigar_pairs = self.parse_cigar(result)
        updated = False
        e.remove_edit(self.consensus_sequence[record_id])

        if self.is_improved_cigar_prefix(cigar_pairs, new_cigar_pairs):
            logging.debug("Keep frame shift as improvement")
            updated = True
            return updated_coordinate_difference, new_cigar_pairs, updated, e
        else:
            logging.debug("Reject frame shift")
            return coordinate_difference, cigar_pairs, updated, e

    def choose_best_frame_shift(self, feature_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs,
                                stop_codons, max_mismatch, coordinate_difference=0):
        """
        Compares the frame shifts obtained by inserting or deleting 1 or 2 letters in the nucleotide query sequence to
        see which if any returns the greatest improvement to the alignment cigar
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record_id:
        :param ref_sequence: amino acid sequence of reference in features
        :param cigar_pairs: best amino acid alignment of query before frame shift
        :param coordinate_difference: amount to offset when applying frame shift (given that we want raw coordinates in
        query consensus sequence and it may have been updated with previously occurring frame shifts
        :return: updated coordinate_difference, updated cigar_pairs, whether updated
        """

        shifts = [("","N"), ("N",""), ("","NN"),("NN","")]
        frame_shift_results = []
        shift_position = self.get_position_for_frame_shift(found_coordinates, record_id, cigar_pairs, stop_codons, max_mismatch)
        logging.debug("Try a frame shift at position %d" % shift_position)
        for shift_from, shift_to in shifts:
            result = self.frame_shift(feature_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs,
                                      shift_from, shift_to, shift_position, coordinate_difference)
            if result[2]:
                frame_shift_results.append(result)

        if len(frame_shift_results) == 0:
            return coordinate_difference, cigar_pairs, False

        logging.debug("Choose winning shift")
        best = 0
        for i,result in enumerate(frame_shift_results):
            if self.is_improved_cigar_prefix(frame_shift_results[best][1], result[1]) \
                    and self.is_longer_cigar_prefix(frame_shift_results[best][1], result[1], max_mismatch):
                best = i
                logging.debug("Override best with %d" %i)

        updated_coordinate_difference, cigar_pairs, updated, edit = frame_shift_results[best]
        logging.debug("Found updated coordinate difference %d" %updated_coordinate_difference)
        edit.apply_edit(self.consensus_sequence[record_id], coordinate_difference)
        self.edits.add_edit(edit)

        return updated_coordinate_difference, cigar_pairs, updated

    def remove_all_edits(self, record_id=0):
        if len(self.edits.edits) == 0:
            return
        record = self.consensus_sequence[record_id]
        self.edits.edits.reverse()
        for edit in self.edits.edits:
            if edit.edit_applied:
                edit.remove_edit(record)

    def str_coordinates(self, coordinates):
        if coordinates is None:
            return coordinates
        else:
            return ",".join([str(i) for i in coordinates])

    def get_in_frame_query_alignment(self, ref_sequence, found_coordinates, record_id, max_mismatch):
        shift = 0
        while shift < 3:
            query_coordinates = (found_coordinates[0]-shift, found_coordinates[1]-shift)
            query_sequence, query_coordinates = self.get_query_sequence(record_id, coordinates=query_coordinates)
            result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
            min_run_length = max_mismatch
            n_runs = self.find_n_runs(query_sequence, min_run_length)
            cigar_pairs = self.parse_cigar(result)
            cigar_length = self.cigar_length(cigar_pairs, max_mismatch, n_runs)
            if cigar_length > min_run_length:
                return query_sequence, query_coordinates, n_runs, cigar_pairs
            shift += 1

        query_sequence, found_coordinates = self.get_query_sequence(record_id, coordinates=found_coordinates)
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        n_runs = self.find_n_runs(query_sequence, min_run_length)
        cigar_pairs = self.parse_cigar(result)

        return query_sequence, found_coordinates, n_runs, cigar_pairs
    
    def discover_frame_shift_edits(self, feature_coordinates, found_coordinates, stop_codons, max_mismatch, record_id=0):
        """
        Gradually introduce frame shifts which improve the amino acid alignment prefix between reference and query
        sequences in an interval
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record_id:
        :return: pairwise alignment with frame shifts added
        """
        logging.debug("Had ref feature coordinates [%s]" % self.str_coordinates(feature_coordinates))
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, shift_into_frame=True)
        logging.debug("Updated ref feature coordinates [%s]" % self.str_coordinates(feature_coordinates))
        logging.debug("Had query feature coordinates [%s]" % self.str_coordinates(found_coordinates))
        query_sequence, found_coordinates, n_runs, cigar_pairs = self.get_in_frame_query_alignment(ref_sequence,
                                                                                                   found_coordinates,
                                                                                                   record_id,
                                                                                                   max_mismatch)
        logging.debug("Updated query feature coordinates [%s]" % self.str_coordinates(found_coordinates))

        coordinate_difference = 0
        while self.cigar_length(cigar_pairs, max_mismatch, n_runs) < len(ref_sequence):
            logging.debug("Cigar shorter than ref: try a frame shift")
            coordinate_difference, cigar_pairs, updated = self.choose_best_frame_shift(feature_coordinates,
                                                                                       found_coordinates,
                                                                                       record_id,
                                                                                       ref_sequence,
                                                                                       cigar_pairs,
                                                                                       stop_codons,
                                                                                       max_mismatch,
                                                                                       coordinate_difference)
            logging.debug("new coordinate difference is %d" %coordinate_difference)
            if not updated:
                break
        logging.debug("Edit list is now: %s" %self.edits)
        self.remove_all_edits(record_id)
        return coordinate_difference

    def add_key_to_coordinate_dict(self, key):
        if key in self.coordinates:
                return
        else:
                self.coordinates[key] = {}

    def save_found_coordinates(self, filepath, write_format='a'):
        with open(filepath,write_format) as f:
            j = json.dumps(self.coordinates)
            f.write(j)
            f.write('\n')

    def run(self, reference_info_filepath, consensus_sequence_filepath, edit_filepath="", stop_codons=["*"],
            max_mismatch=3):
        self.load_input_files(reference_info_filepath, consensus_sequence_filepath, edit_filepath)
        logging.info("Found features: %s " %self.reference_info["references"][self.closest_accession]["locations"])

        for record_id in range(len(self.consensus_sequence)):
            logging.info("Consider consensus sequence %d: %s" %(record_id, self.consensus_sequence[record_id].id))
            for key, value in self.reference_info["references"][self.closest_accession]["locations"].items():
                logging.info("Find edits for %s, %s" %(key,value))
                self.add_key_to_coordinate_dict(key)
                coordinates = (value["start"], self.codon_aware_update(value["start"], value["end"]))
                query_start, query_end = self.identify_feature_coordinates(feature_coordinates=coordinates,
                                                                           record_id=record_id)
                logging.debug("Found feature coordinates [%s]" % self.str_coordinates([query_start, query_end]))
                if not query_end:
                    logging.debug("No good alignment to features coordinates - skip this feature/consensus combination")
                    continue
                logging.debug("Identified features coordinates [%s]" % self.str_coordinates([query_start, query_end]))
                coordinate_difference = self.discover_frame_shift_edits(coordinates, (query_start, query_end),
                                                                        stop_codons, max_mismatch, record_id=record_id)
                logging.info("Total number of discovered edits is %d" %len(self.edits.edits))
                self.coordinates[key][self.consensus_sequence[record_id].id] = {'start': query_start, 'end':
                    query_end + coordinate_difference}



        self.save_found_coordinates(consensus_sequence_filepath + ".coordinates", write_format='w')
        self.edits.save(consensus_sequence_filepath + ".edits", filter_by_applied=False)
