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

class Annotator():
    def __init__(self, closest_accession=None):
        self.consensus_sequence = None
        self.reference_info = None
        self.edits = None
        self.closest_accession = closest_accession
        
    def load_reference_info(self, filepath):
        """
        Load the reference JSON, checking it contains the appropriate fields
        :param filepath:
        :return: dict (from JSON)
        """
        logging.debug("Loading reference JSON %s" % filepath)
        if not os.path.exists(filepath):
            logging.error("Reference filepath %s does not exist!" %filepath)
        with open(filepath) as json_file:
            data = json.load(json_file)
        logging.debug("Checking that JSON has correct format and contains the appropriate fields (sequence and orf) "
                      "for accession %s" %self.closest_accession)
        assert('references' in data.keys())
        assert(self.closest_accession in data['references'].keys())
        assert('sequence' in data['references'][self.closest_accession].keys())
        assert('orf' in data['references'][self.closest_accession].keys())
        return data
        
    def load_consensus_sequence(self, filepath, filetype="fasta"):
        """
        Load consensus FASTA (or other file type) and check there is at least one record
        :param filepath:
        :param filetype: Format of consensus file (accepted by Biopython) default: FASTA
        :return: records
        """
        logging.debug("Loading consensus %s %s" %(filetype,filepath))
        if not os.path.exists(filepath):
            logging.error("Consensus filepath %s does not exist!" %filepath)
        records = list(SeqIO.parse(filepath, filetype))
        logging.debug("The consensus file contains %d records" %len(records))
        assert(len(records) > 0)
        return records
    
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
        
    def get_sequence(self, seq, coordinates=None, offset=None, amino_acid=True):
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
            seq = seq[start:end]
        if offset:
            seq = seq[offset:]
        if amino_acid:
            if len(seq) % 3 == 1:
                seq = seq + "NN"
            elif len(seq) % 3 == 2:
                seq = seq + "NN"
            seq = seq.translate()
        return str(seq)
    
    def get_reference_sequence(self, coordinates=None, offset=None, amino_acid=True, accession=None):
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
        return self.get_sequence(seq, coordinates, offset, amino_acid)
    
    def get_query_sequence(self, record_id=0, coordinates=None, offset=None, amino_acid=True):
        """
        Get the string associated to the given record (within a coordinate range and translated as required)
        :param record_id:
        :param coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param offset:
        :param amino_acid:
        :return:
        """
        seq = self.consensus_sequence[record_id].seq
        return self.get_sequence(seq, coordinates, offset, amino_acid)
    
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
        result = parasail.sw_trace(query_seq, ref_seq, gap_open, gap_extension, matrix)
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
    
    def cigar_length(self,pairs):
        """
        Find the length of aligned sequence until the first insertion/deletion/padding
        :param pairs:
        :return: number
        """
        #pairs = self.parse_cigar(result)
        total = 0
        for c,i in pairs:
            if c in ["M","=","X"]:
                total += i
            elif c in ["I","D","N","S","H","P"]:
                break
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

    def is_improved_cigar_prefix(self, old_cigar_pairs, new_cigar_pairs):
        logging.debug("Old cigar pairs %s" %old_cigar_pairs)
        logging.debug("New cigar pairs %s" %new_cigar_pairs)

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

        #def choose_closest_reference(self):
    #    closest_accession = None
    #    return closest_accession
    
    def identify_orf_coordinates(self, orf_coordinates, record_id=0):
        """
        Find the region in query/consensus sequence which aligns to sequence in reference nucleotide sequence with
        coordinates orf_coordinates
        :param orf_coordinates: 0-based (start,end) in nucleotide sequence (end not included)
        :param record_id: Default 0
        :return:
        """
        ref_sequence = self.get_reference_sequence(orf_coordinates, amino_acid=False)
        query_sequence = self.get_query_sequence(record_id, amino_acid=False)
        result = self.pairwise_ssw_align(ref_sequence, query_sequence)
        return result.read_begin1, result.read_end1+1

    def frame_shift(self, orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs, shift_from,
                    shift_to):
        logging.debug("Frame_shift from '%s' to '%s'" %(shift_from, shift_to))
        cigar_length = self.cigar_length(cigar_pairs)
        e = Edit(record_id, found_coordinates[0] + 3 * (cigar_length), shift_from, shift_to, self.closest_accession,
                 orf_coordinates[0] + 3 * (cigar_length))
        e.apply_edit(self.consensus_sequence[record_id])
        coordinate_difference = len(shift_to) - len(shift_from)
        updated_found_coordinates = [found_coordinates[0], found_coordinates[1] + coordinate_difference]
        query_sequence = self.get_query_sequence(record_id, coordinates=updated_found_coordinates)
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        new_cigar_pairs = self.parse_cigar(result)
        updated = False
        e.remove_edit(self.consensus_sequence[record_id])

        if self.is_improved_cigar_prefix(cigar_pairs, new_cigar_pairs):
            logging.debug("Keep frame shift as improvement")
            updated = True
            return updated_found_coordinates, new_cigar_pairs, updated, e
        else:
            logging.debug("Reject frame shift")
            return found_coordinates, cigar_pairs, updated, e

    def choose_best_frame_shift(self, orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs):

        shifts = [("","N"), ("N",""), ("","NN"),("NN","")]
        frame_shift_results = []
        for shift_from, shift_to in shifts:
            result = self.frame_shift(orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_pairs,
                                      shift_from, shift_to)
            if result[2]:
                frame_shift_results.append(result)

        if len(frame_shift_results) == 0:
            return found_coordinates, cigar_pairs, False

        logging.debug("Choose winning shift")
        best = 0
        for i,result in enumerate(frame_shift_results):
            if self.is_improved_cigar_prefix(frame_shift_results[best][1], result[1]):
                best = i
                logging.debug("Override best with %d" %i)

        found_coordinates, cigar_pairs, updated, edit = frame_shift_results[best]
        edit.apply_edit(self.consensus_sequence[record_id])
        self.edits.add_edit(edit)

        return found_coordinates, cigar_pairs, updated
    
    def discover_edits(self, orf_coordinates, found_coordinates, record_id=0):
        ref_sequence = self.get_reference_sequence(orf_coordinates)
        query_sequence = self.get_query_sequence(record_id, coordinates=found_coordinates)
        
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        cigar_pairs = self.parse_cigar(result)
        while self.cigar_length(cigar_pairs) < len(ref_sequence):
            logging.debug("Cigar shorter than ref: try a frame shift")
            found_coordinates, cigar_pairs, updated = self.choose_best_frame_shift(orf_coordinates, found_coordinates,
                                                                                   record_id, ref_sequence, cigar_pairs)
            if not updated:
                break
        logging.debug("Edit list is now: %s" %self.edits)
        return result

    def run(self, reference_info_filepath, consensus_sequence_filepath, edit_filepath=""):
        self.load_input_files(reference_info_filepath, consensus_sequence_filepath, edit_filepath)
        print(self.reference_info["references"][self.closest_accession]["orf"])
        for key, value in self.reference_info["references"][self.closest_accession]["orf"].items():
            print(key,value)
            coordinates = (value["start"], value["end"])
            query_start, query_end = self.identify_orf_coordinates(orf_coordinates=coordinates)
            print(query_start, query_end)
            result = self.discover_edits(coordinates, (query_start, query_end))
