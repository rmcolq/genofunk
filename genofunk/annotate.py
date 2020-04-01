import os
import sys
import logging
from Bio import SeqIO
import json  
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import parasail

from genofunk import editfile
from genofunk.parasail_utils import *
from genofunk.sequence_utils import *

EditFile = editfile.EditFile
Edit = editfile.Edit


class Annotate:
    def __init__(self, closest_accession=None):
        self.consensus_sequence = None
        self.reference_info = None
        self.edits = None
        self.coordinates = {}
        self.closest_accession = closest_accession
        self.problematic = set()
        
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
        self.consensus_sequence = SeqIO.index(filepath, filetype)
        logging.debug("The consensus file contains %d records" %len(self.consensus_sequence))
        assert(len(self.consensus_sequence) > 0)

    def apply_loaded_edits(self):
        """
        Apply the edits to the consensus nucleotide sequences in place (in reverse order to avoid offset errors)
        :return:
        """
        if len(self.edits.edits) == 0:
            return
        self.edits.sort(reverse=True, seq_position=True)
        for edit in self.edits.edits:
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
        self.load_consensus_sequence(consensus_filepath)
        self.reference_info = self.load_reference_info(reference_filepath)
        self.edits = EditFile(edit_filepath)
        self.apply_loaded_edits()
    
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
        return get_sequence(seq, coordinates, shift_into_frame, offset, amino_acid)

    def get_query_sequence(self, record, coordinates=None, shift_into_frame=False, offset=None, amino_acid=True):
        """
        Get the string associated to the given record (within a coordinate range and translated as required)
        :param record or record_id:
        :param coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param offset:
        :param amino_acid:
        :return:
        """
        if type(record) != SeqRecord:
            seq = self.consensus_sequence[record].seq
            return get_sequence(seq, coordinates, shift_into_frame, offset, amino_acid)
        return get_sequence(record.seq, coordinates, shift_into_frame, offset, amino_acid)

    def get_query_length(self, record_id=0):
        """
        Get the length of the string associated to the given record
        :param record_id:
        :return: length
        """
        seq = self.consensus_sequence[record_id].seq
        return len(seq)

    def update_coordinates_if_n_runs(self, query_start, query_end, query_sequence, ref_match_start, ref_match_end,
                                     feature_length, result):
        first_run_n = find_run_n(query_sequence)
        last_run_n = rfind_run_n(query_sequence)
        logging.debug("Found first run Ns %i,%i and last run Ns %i,%i" % (first_run_n[0], first_run_n[1], last_run_n[0],
                                                                          last_run_n[1]))

        ref_match_length = ref_match_end - ref_match_start
        tolerance = 32

        n_run_closer_to_start_than_end = first_run_n[0] < query_end - first_run_n[1]
        strong_start_match = False
        strong_end_match = False
        cigar_pairs = parse_cigar_pairs(result)
        if len(cigar_pairs) > 2 and cigar_pairs[1][0] == '=' and cigar_pairs[1][1] > tolerance:
            strong_start_match = True
        elif len(cigar_pairs) > 2 and cigar_pairs[-1][0] == '=' and cigar_pairs[-1][1] > tolerance:
            strong_end_match = True

        if first_run_n[0] == -1 and last_run_n[0] == -1:
            logging.debug("No run Ns - did not update coordinates")
        elif first_run_n[0] == last_run_n[0] \
                and (strong_end_match or (not strong_start_match and n_run_closer_to_start_than_end)) \
                and query_start + first_run_n[0] <= query_end - feature_length + tolerance:
            #query_start = codon_aware_update(query_start, query_end - feature_length)
            query_start = query_end - feature_length
            logging.debug("Only one run Ns at start")
        elif first_run_n[0] == last_run_n[0] \
                and (strong_start_match or (not strong_end_match and not n_run_closer_to_start_than_end)) \
                and query_start + first_run_n[1] >= query_start + feature_length - tolerance:
            #query_end = codon_aware_update(query_start, query_start + feature_length)
            query_end = query_start + feature_length
            logging.debug("Only one run Ns at end")
        elif first_run_n[0] != last_run_n[0] \
                and last_run_n[0] - first_run_n[1] > feature_length:
            #query_end = codon_aware_update(query_start, query_start + last_run_n[1])
            query_end = query_start + last_run_n[0] - 1
            #query_start = query_start + first_run_n[0]
            query_start = codon_aware_update(query_start, query_start + first_run_n[0])
            logging.debug("Two runs Ns at start and end with long gap between")
        elif first_run_n[0] != last_run_n[0] \
                and last_run_n[1] - first_run_n[0] > feature_length:
            #query_end = codon_aware_update(query_start, query_start + last_run_n[1])
            query_end = query_start + last_run_n[1]
            #query_start = query_start + first_run_n[0]
            query_start = codon_aware_update(query_start, query_start + first_run_n[0])
            logging.debug("Two runs Ns at start and end")
        elif first_run_n[0] != last_run_n[0] \
                and (strong_end_match or first_run_n[0] <= query_end - feature_length <= first_run_n[1]):
            #query_start = codon_aware_update(query_start, query_end - feature_length)
            query_start = query_end - feature_length
            logging.debug("Two runs Ns, including at start")
        elif first_run_n[0] != last_run_n[0] \
                and (strong_start_match or last_run_n[0] <= feature_length <= last_run_n[1]):
            #query_end = codon_aware_update(query_start, query_start + feature_length)
            query_end = query_start + feature_length
            logging.debug("Two runs Ns, including at end")
        else:
            logging.debug("Did not update coordinates")
        return query_start, query_end

    def identify_feature_coordinates(self, feature_coordinates, record_id=0, min_score = 300,
                                     max_query_start_offset = 200, max_nucleotide_length_difference = 10):
        """
        Find the region in query/consensus sequence which aligns to sequence in reference nucleotide sequence with
        coordinates feature_coordinates
        :param feature_coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param record_id: Default 0
        :return:
        """
        logging.debug("Have ref feature coordinates %s" % str_coordinates(feature_coordinates))
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, amino_acid=False)
        if len(ref_sequence) < 60:
            min_score = 5 * len(ref_sequence)

        query_seq_length = self.get_query_length(record_id)
        offset = int((feature_coordinates[0] - max_query_start_offset)*0.99)
        if offset > query_seq_length:
            offset = 0
        query_end = codon_aware_update(offset, int((feature_coordinates[1] + max_query_start_offset)*1.01))
        query_sequence, coordinates = self.get_query_sequence(record_id,
                                                              coordinates=(offset, query_end),
                                                              amino_acid=False)
        logging.debug("Start with query feature coordinates %s" % str_coordinates(coordinates))
        result = pairwise_sw_align(ref_sequence, query_sequence)
        if result and result.score > min_score and cigar_has_min_matches(result, min_matches=10):
            (ref_begin, ref_end, read_begin, read_end) = get_alignment_start_end(result)
            logging.debug("Found score %s and cigar %s" % (result.score, result.cigar.decode))
            logging.debug("Found query match start and end %s, %s" % (read_begin, read_end))
            logging.debug("Found ref match start and end %s, %s" % (ref_begin, ref_end))

            query_start = offset + read_begin - ref_begin
            if query_start > query_seq_length:
                query_start = 0
            logging.debug("Update query start to %i" % query_start)
            feature_length = feature_coordinates[1] - feature_coordinates[0]
            #query_end = min(query_seq_length, codon_aware_update(query_start, offset + read_end + feature_length - ref_end))
            query_end = min(query_seq_length, offset + read_end + feature_length - ref_end)
            logging.debug("Update query end to %i" % query_end)
            found_length = query_end - query_start
            query_sequence, coordinates = self.get_query_sequence(record_id,
                                                                  coordinates=(query_start, query_end),
                                                                  amino_acid=False)
            logging.debug("Found length %i and feature length %i" %(found_length, feature_length))

            if found_length - feature_length > 0:
                updated = True
                attempt = 0
                while found_length - feature_length > 0 and updated and attempt < 3:
                    query_start, query_end = self.update_coordinates_if_n_runs(query_start, query_end, query_sequence,
                                                                               ref_begin, ref_end, feature_length, result)
                    new_found_length = query_end - query_start
                    if new_found_length == found_length:
                        updated = False
                    else:
                        query_sequence, coordinates = self.get_query_sequence(record_id,
                                                                              coordinates=(query_start, query_end),
                                                                              amino_acid=False)
                        query_start, query_end = coordinates
                    found_length = new_found_length
                    attempt += 1
                    logging.debug("At attempt %i, coordinates updated was %s and new found length is now %i"
                                  %(attempt, str(updated), found_length))
                query_sequence, coordinates = self.get_query_sequence(record_id,
                                                                      coordinates=(query_start, query_end),
                                                                      #shift_into_frame=True,
                                                                      amino_acid=False)
                query_start, query_end = coordinates
                logging.debug("Update query start, end to %i, %i" % (query_start, query_end))
                logging.debug("Found length %i and feature length %i" % (found_length, feature_length))

            return query_start, query_end
        else:
            if result:
                logging.debug("Found result, but score was %i (not > %i) or has insufficient match bases" % (result.score, min_score))
            else:
                logging.debug("No result")
                logging.debug("Ref sequence %s" % ref_sequence)
                logging.debug("Qry sequence %s" % query_sequence)

            return None, None

    def get_n_codon_positions(self, record, found_coordinates, min_frame_shift_position):
        positions = []
        query_sequence, coordinates = self.get_query_sequence(record, coordinates=found_coordinates, amino_acid=False)

        next = "N"
        min_pos = 3*min_frame_shift_position
        attempts = 0
        while next == "N" and attempts < 100:
            position = query_sequence.find("NNN", min_pos)
            if 0 < position < len(query_sequence) - 1:
                while position < len(query_sequence) - 1 and query_sequence[position + 1] == "N":
                    position += 1
                if position + 1 > len(query_sequence) - 1:
                    position = -1
                    break
                next = query_sequence[position + 1]
            elif position > 0 and position % 3 != 0:
                next = "N"
            else:
                break
            min_pos = position + 1
            logging.debug("%i %s %i" %(position, next, min_pos))
            attempts += 1
        if position > 0 and position % 3 == 0:
            positions.append(int(position/3))
            logging.debug("Found NNN position %i" %int(position/3))
        return positions

    def get_cigar_based_positions(self, query_sequence, cigar_pairs, max_mismatch, min_frame_shift_position):
        positions = []
        min_run_length = max_mismatch
        n_runs = find_n_runs(query_sequence, min_run_length)
        position = get_position_first_indel_or_mismatch_in_cigar(cigar_pairs, max_mismatch, n_runs, min_position=min_frame_shift_position)
        if position > min_frame_shift_position:
            positions.append(position)
        return positions

    def get_position_for_frame_shift(self, found_coordinates, record, cigar_pairs, stop_codons, max_mismatch,
                                     min_frame_shift_position=0):
        #logging.debug("Get position for frame shift")
        positions = []
        query_sequence, coordinates = self.get_query_sequence(record, coordinates=found_coordinates)
        #logging.debug("Have query sequence %s and coordinates %s" % (query_sequence, coordinates))

        positions.extend(get_stop_codon_positions(query_sequence, stop_codons, min_frame_shift_position))
        #logging.debug(positions)
        positions.extend(self.get_n_codon_positions(record, found_coordinates, min_frame_shift_position))
        #logging.debug(positions)
        positions.extend(self.get_cigar_based_positions(query_sequence, cigar_pairs, max_mismatch, min_frame_shift_position))
        #logging.debug(positions)

        if len(positions) == 0:
            return None
        return min(positions)


    def frame_shift(self, feature_coordinates, found_coordinates, record, ref_sequence, cigar_pairs, stop_codons,
                    shift_from, shift_to, shift_position, include_compensatory, coordinate_difference=0):
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

        record_name = record.id
        old_query_sequence, coordinates = self.get_query_sequence(record, coordinates=found_coordinates)
        e = Edit(record_name, found_coordinates[0] + 3 * (shift_position), shift_from, shift_to, self.closest_accession,
                 feature_coordinates[0] + 3 * (shift_position))
        record, applied = e.apply_edit(record, coordinate_difference)
        updated_coordinate_difference = coordinate_difference
        if applied:
            updated_coordinate_difference += len(shift_to) - len(shift_from)
        updated_found_coordinates = [found_coordinates[0], found_coordinates[1] + updated_coordinate_difference]
        query_sequence, coordinates = self.get_query_sequence(record, coordinates=updated_found_coordinates)
        logging.debug("Found query sequence %s and coordinates %s" %(query_sequence, coordinates))
        logging.debug("Have ref sequence %s" % ref_sequence)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        new_cigar_pairs = parse_cigar_pairs(result)
        updated = False
        record = e.remove_edit(record)

        logging.debug("include compensatory %s, has_fewer_stop_codons %s, is improved_cigar %s"
                      %(include_compensatory, has_fewer_stop_codons(old_query_sequence, query_sequence, stop_codons),
                        is_improved_cigar(cigar_pairs, new_cigar_pairs)))
        if (include_compensatory or has_fewer_stop_codons(old_query_sequence, query_sequence, stop_codons)) and \
                is_improved_cigar(cigar_pairs, new_cigar_pairs):
            logging.debug("Keep frame shift as improvement")
            updated = True
            return updated_coordinate_difference, new_cigar_pairs, updated, e
        else:
            logging.debug("Reject frame shift")
            return coordinate_difference, cigar_pairs, updated, e

    def choose_best_frame_shift(self, feature_coordinates, found_coordinates, record, ref_sequence, cigar_pairs,
                                stop_codons, max_mismatch, include_compensatory, coordinate_difference=0,
                                shift_position=0):
        """
        Compares the frame shifts obtained by inserting or deleting 1 or 2 letters in the nucleotide query sequence to
        see which if any returns the greatest improvement to the alignment cigar
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record:
        :param ref_sequence: amino acid sequence of reference in features
        :param cigar_pairs: best amino acid alignment of query before frame shift
        :param coordinate_difference: amount to offset when applying frame shift (given that we want raw coordinates in
        query consensus sequence and it may have been updated with previously occurring frame shifts
        :return: updated coordinate_difference, updated cigar_pairs, whether updated
        """

        shifts = [("","N"), ("N",""), ("","NN"),("NN","")]
        frame_shift_results = []

        for shift_from, shift_to in shifts:
            result = self.frame_shift(feature_coordinates, found_coordinates, record, ref_sequence, cigar_pairs,
                                      stop_codons, shift_from, shift_to, shift_position, include_compensatory,
                                      coordinate_difference)
            if result[2]:
                frame_shift_results.append(result)

        if len(frame_shift_results) == 0:
            return coordinate_difference, cigar_pairs, False, record

        logging.debug("Choose winning shift")
        best = 0
        for i,result in enumerate(frame_shift_results):
            if is_improved_cigar(frame_shift_results[best][1], result[1], consider_if_frameshift_added=False) \
                    and is_longer_cigar_prefix(frame_shift_results[best][1], result[1], max_mismatch):
                best = i
                logging.debug("Override best with %d" %i)

        updated_coordinate_difference, cigar_pairs, updated, edit = frame_shift_results[best]
        logging.debug("Found updated coordinate difference %d" %updated_coordinate_difference)
        record, applied = edit.apply_edit(record, coordinate_difference)
        self.edits.add_edit(edit)

        return updated_coordinate_difference, cigar_pairs, updated, record

    def remove_all_edits(self, record):
        if len(self.edits.edits) == 0:
            return
        self.edits.sort(reverse=True, seq_position=True)
        for edit in self.edits.edits:
            if edit.edit_applied:
                record = edit.remove_edit(record)

    def get_in_frame_query_alignment(self, ref_sequence, found_coordinates, record_id, max_mismatch):
        shift = 0
        #logging.debug("start with found coordinates %s" % str_coordinates(found_coordinates))

        while shift < 3:
            query_coordinates = (found_coordinates[0]-shift, found_coordinates[1]-shift)
       #     logging.debug("query coordinates %s" % str_coordinates(query_coordinates))
            query_sequence, query_coordinates = self.get_query_sequence(record_id, coordinates=query_coordinates)
       #     logging.debug("getting in frame query alignment with query_sequence %s and ref sequence %s" %(query_sequence, ref_sequence))
            result = pairwise_nw_trace_align(ref_sequence, query_sequence)
            min_run_length = max_mismatch
            n_runs = find_n_runs(query_sequence, min_run_length)
            cigar_pairs = parse_cigar_pairs(result)
            cigar_length = get_position_first_indel_or_mismatch_in_cigar(cigar_pairs, max_mismatch, n_runs)
            if cigar_length > min_run_length:
                return query_sequence, query_coordinates, n_runs, cigar_pairs
            shift += 1

        query_sequence, found_coordinates = self.get_query_sequence(record_id, coordinates=found_coordinates)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        n_runs = find_n_runs(query_sequence, min_run_length)
        cigar_pairs = parse_cigar_pairs(result)

        return query_sequence, found_coordinates, n_runs, cigar_pairs
    
    def discover_frame_shift_edits(self, feature_coordinates, found_coordinates, stop_codons, max_mismatch,
                                   include_compensatory, record_id=0):
        """
        Gradually introduce frame shifts which improve the amino acid alignment prefix between reference and query
        sequences in an interval
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record_id:
        :return: pairwise alignment with frame shifts added
        """
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates)
        query_sequence, found_coordinates, n_runs, cigar_pairs = self.get_in_frame_query_alignment(ref_sequence,
                                                                                                   found_coordinates,
                                                                                                   record_id,
                                                                                                   max_mismatch)
        logging.debug("found coordinates %s and query_sequence %s" %(str_coordinates(found_coordinates),query_sequence))
        coordinate_difference = 0
        record = self.consensus_sequence[record_id]
        min_frame_shift_position = 0
        candidate_shift_position = self.get_position_for_frame_shift(found_coordinates, record, cigar_pairs,
                                                                    stop_codons, max_mismatch, min_frame_shift_position)
        attempt = 0
        while candidate_shift_position is not None and candidate_shift_position < len(ref_sequence) and attempt < 100:
            logging.debug("Cigar shorter than ref: try a frame shift")
            coordinate_difference, cigar_pairs, updated, record = \
                                                          self.choose_best_frame_shift(feature_coordinates,
                                                                                       found_coordinates,
                                                                                       record,
                                                                                       ref_sequence,
                                                                                       cigar_pairs,
                                                                                       stop_codons,
                                                                                       max_mismatch,
                                                                                       include_compensatory,
                                                                                       coordinate_difference,
                                                                                       candidate_shift_position)
            #logging.debug("new coordinate difference is %d" %coordinate_difference)
            #logging.debug("Checked found coordinates %s" % str_coordinates(found_coordinates))
            min_frame_shift_position = candidate_shift_position
            candidate_shift_position = self.get_position_for_frame_shift(found_coordinates, record, cigar_pairs,
                                                                    stop_codons, max_mismatch, min_frame_shift_position)
            attempt += 1
            if not updated:
                break
        logging.debug("Edit list is now: %s" %self.edits)
        self.remove_all_edits(record)
        return coordinate_difference, found_coordinates

    def discover_sequence_edits(self, feature_coordinates, found_coordinates, stop_codons, max_mismatch, record_id=0):
        """
        Find positions where there is a sequence mismatch and store them, assuming reference is correct. When combined
        with information during the merge step, will result in a list of edits which are unique to this sequence.
        :param feature_coordinates: coordinates of features in reference sequence
        :param found_coordinates: corresponding coordinates for features in query sequence
        :param record_id:
        :return: pairwise alignment with frame shifts added
        """
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, amino_acid=False)
        query_sequence, found_coordinates = self.get_query_sequence(found_coordinates, amino_acid=False)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = parse_cigar_pairs(result)

        record = self.consensus_sequence[record_id]
        min_edit_position = 0
        candidate_edit_position, ref_position, candidate_length = get_position_next_short_mismatch_in_cigar(pairs,
                                                                                                    max_mismatch,
                                                                                                    min_edit_position)
        attempt = 0
        while candidate_edit_position < len(ref_sequence) and attempt < 100:
            logging.debug("Found edit at position %i" %candidate_edit_position)
            shift_from = query_sequence[candidate_edit_position:candidate_edit_position + candidate_length + 1]
            shift_to = ref_sequence[ref_position:ref_position + candidate_length + 1]
            e = Edit(record_id, found_coordinates[0] + candidate_edit_position, shift_from, shift_to,
                     self.closest_accession,
                     feature_coordinates[0] + ref_position)
            if "N" not in shift_from:
                self.edits.add_edit(e)
            min_edit_position = candidate_edit_position
            candidate_edit_position, ref_position, candidate_length = get_position_next_short_mismatch_in_cigar(pairs,
                                                                                                    max_mismatch,
                                                                                                    min_edit_position)
            attempt += 1

        logging.debug("Edit list is now: %s" %self.edits)

    def save_found_coordinates(self, filepath, write_format='a'):
        with open(filepath,write_format) as f:
            j = json.dumps(self.coordinates)
            f.write(j)
            f.write('\n')

    def check_have_open_reading_frame_ref(self, ref_coordinate_pairs):
        all_ref_coordinates = []
        for pair in ref_coordinate_pairs:
            all_ref_coordinates.extend(pair)

        for i in range(1):
            all_ref_coordinates[-1 ] += i
            ref_sequence, coordinates = self.get_reference_sequence(coordinates=all_ref_coordinates)
            if is_open_reading_frame(ref_sequence, allow_missing=False):
                logging.debug("Reference sequence %s with coordinates %s has GOOD open reading frame"
                              % (ref_sequence, str_coordinates(all_ref_coordinates)))
                ref_coordinate_pairs[-1][-1] += i
                return ref_coordinate_pairs
        logging.debug("Reference sequence %s with coordinates %s has bad open reading frame"
                          %(ref_sequence, str_coordinates(all_ref_coordinates)))
        sys.exit()

    def check_have_open_reading_frame_query(self, record_id, query_coordinate_pairs, allow_stop_codons_in_middle):
        all_query_coordinates = []
        for pair in query_coordinate_pairs:
            all_query_coordinates.append(pair['start'])
            all_query_coordinates.append(pair['end'])
        record = self.consensus_sequence[record_id]
        record,coordinate_difference = apply_edits_in_range(self.edits, record, coordinates=all_query_coordinates)

        query_sequence, coordinates = self.get_query_sequence(record, coordinates=all_query_coordinates)

        if not is_open_reading_frame(query_sequence, allow_stop_codons_in_middle=allow_stop_codons_in_middle):
            logging.debug("After finding edits, still have bad open reading frame for query sequence %s with "
                          "coordinates %s" %(query_sequence, str_coordinates(all_query_coordinates)))
            self.problematic.add(record_id)
        else:
            logging.debug("Good reading frame for query sequence %s with "
                          "coordinates %s" %(query_sequence, str_coordinates(all_query_coordinates)))
        self.remove_all_edits(record)

    def run(self, reference_info_filepath, consensus_sequence_filepath, edit_filepath="", stop_codons=["*"],
            max_mismatch=3, include_compensatory=False, min_seq_length=28000, discover_frame_shifts=False,
            allow_stop_codons_in_middle=True):
        self.load_input_files(reference_info_filepath, consensus_sequence_filepath, edit_filepath)
        logging.info("Found features: %s " %self.reference_info["references"][self.closest_accession]["locations"])

        for record_id in self.consensus_sequence:
            logging.info("Consider consensus sequence %s" % record_id)
            if self.get_query_length(record_id) < min_seq_length:
                logging.info("Skip sequence as shorter than min_seq_length %d" % min_seq_length)
                continue
            for key, value in self.reference_info["references"][self.closest_accession]["locations"].items():
                logging.info("Find edits for %s, %s" %(key,value))
                add_key_to_coordinate_dict(self.coordinates, key)

                feature_coordinate_pairs = get_coordinates_from_json(value, pairs=True)
                self.check_have_open_reading_frame_ref(feature_coordinate_pairs)

                query_coordinate_pairs = []
                for coordinates in feature_coordinate_pairs:
                    query_coordinates = self.identify_feature_coordinates(feature_coordinates=coordinates,
                                                                               record_id=record_id)
                    if not query_coordinates[1]:
                        logging.debug("No good alignment to features coordinates - skip this feature/consensus combination")
                        self.problematic.add(record_id)
                        continue
                    logging.debug("Identified features coordinates %s" % str_coordinates(query_coordinates))

                    if discover_frame_shifts:
                        coordinate_difference, query_coordinates = self.discover_frame_shift_edits(coordinates, query_coordinates,
                                                                                stop_codons, max_mismatch, include_compensatory, record_id=record_id)
                        logging.info("Total number of discovered edits is %d" %len(self.edits.edits))
                        logging.debug("Check features coordinates %s" % str_coordinates(query_coordinates))

                    query_coordinate_pairs.append({'start': query_coordinates[0], 'end': query_coordinates[1] + coordinate_difference})

                self.check_have_open_reading_frame_query(record_id, query_coordinate_pairs, allow_stop_codons_in_middle)

                if len(query_coordinate_pairs) == 1:
                    self.coordinates[key][self.consensus_sequence[record_id].id] = query_coordinate_pairs[0]
                else:
                    self.coordinates[key][self.consensus_sequence[record_id].id] = {"join": query_coordinate_pairs}

        self.save_found_coordinates(consensus_sequence_filepath + ".coordinates", write_format='w')
        self.edits.save(consensus_sequence_filepath + ".edits", filter_by_applied=False)
        self.consensus_sequence.close()

        if len(self.problematic) > 0:
            logging.info("A number of records were problematic and after auto-discovering edits they did not have the "
                         "appropriate format for one or more features. This could mean the feature was not found,"
                         "or that the sequence with applied edits did not start with an M and end with a stop, or may"
                         "have had stops in the middle. Please check manually and update editfile.")
            logging.info("These were:")
            for id in self.problematic:
                logging.info(id)
