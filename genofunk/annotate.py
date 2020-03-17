import os
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
            end = min(end, len(seq))
            seq = seq[start:end]
            coordinates = (start, end)
        if offset:
            seq = seq[offset:]
        if shift_into_frame:
            seq, coordinates = shift_nucleotide_sequence_into_frame(seq, coordinates)
        else:
            seq = make_sequence_length_divide_by_3(seq)
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

    def get_query_length(self, record_id=0):
        """
        Get the length of the string associated to the given record
        :param record_id:
        :return: length
        """
        seq = self.consensus_sequence[record_id].seq
        return len(seq)

    def update_coordinates_if_n_runs(self, query_start, query_end, query_sequence, ref_match_start, ref_match_end,
                                     feature_length):
        first_run_n = find_run_n(query_sequence)
        last_run_n = rfind_run_n(query_sequence)
        logging.debug("Found first run Ns %i,%i and last run Ns %i,%i" % (first_run_n[0], first_run_n[1], last_run_n[0],
                                                                          last_run_n[1]))

        ref_match_length = ref_match_end - ref_match_start
        tolerance = 32
        if first_run_n[0] == last_run_n[0] \
                and first_run_n[0] < query_end - first_run_n[1] \
                and query_start + first_run_n[0] <= query_end - feature_length + tolerance:
            query_start = codon_aware_update(query_start, query_end - feature_length)
            logging.debug("Only one run Ns at start")
        elif first_run_n[0] == last_run_n[0] \
                and first_run_n[0] > query_end - first_run_n[1] \
                and query_start + first_run_n[1] >= query_start + feature_length - tolerance:
            query_end = codon_aware_update(query_start, query_start + feature_length)
            logging.debug("Only one run Ns at end")
        elif first_run_n[0] != last_run_n[0] and last_run_n[1] - first_run_n[0] > feature_length:
            query_start = codon_aware_update(query_start, query_start + first_run_n[0])
            query_end = codon_aware_update(query_start, query_start + last_run_n[1])
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
        logging.debug("Have ref feature coordinates %s" % str_coordinates(feature_coordinates))
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, amino_acid=False)

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
            logging.debug("Found query start and end %s, %s" % (read_begin, read_end))
            logging.debug("Found ref start and end %s, %s" % (ref_begin, ref_end))


            query_start = offset + read_begin - ref_begin
            if query_start > query_seq_length:
                query_start = 0
            logging.debug("Update query start to %i" % query_start)
            feature_length = feature_coordinates[1] - feature_coordinates[0]
            query_end = min(query_seq_length, codon_aware_update(query_start, offset + read_end + feature_length - ref_end))
            logging.debug("Update query end to %i" % query_end)
            found_length = query_end - query_start
            query_sequence, coordinates = self.get_query_sequence(record_id,
                                                                  coordinates=(query_start, query_end),
                                                                  amino_acid=False)
            logging.debug("Query sequence %s" %query_sequence)
            logging.debug("Found length %i and feature length %i" %(found_length, feature_length))

            if found_length - feature_length > 0:
                query_start, query_end = self.update_coordinates_if_n_runs(query_start, query_end, query_sequence,
                                                                           ref_begin, ref_end,
                                                                           feature_length)
                logging.debug("Update query start, end to %i, %i" % (query_start, query_end))
                found_length = query_end - query_start
                #query_sequence, coordinates = self.get_query_sequence(record_id, coordinates=(query_start, query_end),
                #                                                      amino_acid=False)
                #query_start, query_end = coordinates
                #logging.debug("Query sequence %s" % query_sequence)
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
        n_runs = find_n_runs(query_sequence,min_run_length)
        cigar_length = get_cigar_length(cigar_pairs, max_mismatch, n_runs)
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
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        new_cigar_pairs = parse_cigar_pairs(result)
        updated = False
        e.remove_edit(self.consensus_sequence[record_id])

        if is_improved_cigar_prefix(cigar_pairs, new_cigar_pairs):
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
            if is_improved_cigar_prefix(frame_shift_results[best][1], result[1], consider_if_frameshift_added=False) \
                    and is_longer_cigar_prefix(frame_shift_results[best][1], result[1], max_mismatch):
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

    def get_in_frame_query_alignment(self, ref_sequence, found_coordinates, record_id, max_mismatch):
        shift = 0
        while shift < 3:
            query_coordinates = (found_coordinates[0]-shift, found_coordinates[1]-shift)
            query_sequence, query_coordinates = self.get_query_sequence(record_id, coordinates=query_coordinates)
            result = pairwise_nw_trace_align(ref_sequence, query_sequence)
            min_run_length = max_mismatch
            n_runs = find_n_runs(query_sequence, min_run_length)
            cigar_pairs = parse_cigar_pairs(result)
            cigar_length = get_cigar_length(cigar_pairs, max_mismatch, n_runs)
            if cigar_length > min_run_length:
                return query_sequence, query_coordinates, n_runs, cigar_pairs
            shift += 1

        query_sequence, found_coordinates = self.get_query_sequence(record_id, coordinates=found_coordinates)
        result = pairwise_nw_trace_align(ref_sequence, query_sequence)
        n_runs = find_n_runs(query_sequence, min_run_length)
        cigar_pairs = parse_cigar_pairs(result)

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
        logging.debug("Had ref feature coordinates [%s]" % str_coordinates(feature_coordinates))
        ref_sequence, feature_coordinates = self.get_reference_sequence(feature_coordinates, shift_into_frame=True)
        logging.debug("Updated ref feature coordinates [%s]" % str_coordinates(feature_coordinates))
        logging.debug("Had query feature coordinates [%s]" % str_coordinates(found_coordinates))
        query_sequence, found_coordinates, n_runs, cigar_pairs = self.get_in_frame_query_alignment(ref_sequence,
                                                                                                   found_coordinates,
                                                                                                   record_id,
                                                                                                   max_mismatch)
        logging.debug("Updated query feature coordinates [%s]" % str_coordinates(found_coordinates))

        coordinate_difference = 0
        while get_cigar_length(cigar_pairs, max_mismatch, n_runs) < len(ref_sequence):
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
            max_mismatch=3, min_seq_length=28000):
        self.load_input_files(reference_info_filepath, consensus_sequence_filepath, edit_filepath)
        logging.info("Found features: %s " %self.reference_info["references"][self.closest_accession]["locations"])

        for record_id in range(len(self.consensus_sequence)):
            logging.info("Consider consensus sequence %d: %s" %(record_id, self.consensus_sequence[record_id].id))
            if self.get_query_length(record_id) < min_seq_length:
                logging.info("Skip sequence as shorter than min_seq_length %d" % min_seq_length)
                continue
            for key, value in self.reference_info["references"][self.closest_accession]["locations"].items():
                logging.info("Find edits for %s, %s" %(key,value))
                self.add_key_to_coordinate_dict(key)
                coordinates = (value["start"], codon_aware_update(value["start"], value["end"]))
                query_start, query_end = self.identify_feature_coordinates(feature_coordinates=coordinates,
                                                                           record_id=record_id)
                #logging.debug("Found feature coordinates [%s]" % str_coordinates([query_start, query_end]))
                if not query_end:
                    logging.debug("No good alignment to features coordinates - skip this feature/consensus combination")
                    continue
                logging.debug("Identified features coordinates [%s]" % str_coordinates([query_start, query_end]))
                coordinate_difference = self.discover_frame_shift_edits(coordinates, (query_start, query_end),
                                                                        stop_codons, max_mismatch, record_id=record_id)
                logging.info("Total number of discovered edits is %d" %len(self.edits.edits))
                self.coordinates[key][self.consensus_sequence[record_id].id] = {'start': query_start, 'end':
                    query_end + coordinate_difference}

        self.save_found_coordinates(consensus_sequence_filepath + ".coordinates", write_format='w')
        self.edits.save(consensus_sequence_filepath + ".edits", filter_by_applied=False)
