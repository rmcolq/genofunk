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
        # tests:
        # what if no file
        # simple case check
        # expected attributes (reference,sequence), (genes,(start,end,strand))
        with open(filepath) as json_file:
            data = json.load(json_file)
        return data
        
    def load_consensus_sequence(self, filepath, filetype="fasta"):
        # tests:
        # what if no file
        # what if file wrong format
        # simple case check
        #record_dict = SeqIO.to_dict(SeqIO.parse(filepath, filetype))
        records = list(SeqIO.parse(filepath, filetype))
        return records
    
    def load_input_files(self, consensus_filepath, reference_filepath, edit_filepath = ""):
        # tests:
        # filepath doesn't exist
        # filepath does exist
        self.consensus_sequence = self.load_consensus_sequence(consensus_filepath)
        self.reference_info = self.load_reference_info(reference_filepath)
        self.edits = EditFile(edit_filepath)
        
    def get_sequence(self, seq, coordinates=None, offset=None, amino_acid=True):
        if coordinates:
            #0 based, not including end coordinate
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
        if not accession:
            accession = self.closest_accession
        seq = Seq(self.reference_info['references'][accession]['sequence'])
        return self.get_sequence(seq, coordinates, offset, amino_acid)
    
    def get_query_sequence(self, record_id=0, coordinates=None, offset=None, amino_acid=True):
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
            
    def pairwise_ssw_align(self, ref_seq, query_seq):
        #result = parasail.sw_trace(query_seq, ref_seq, 10, 1, parasail.blosum62)
        result = parasail.ssw(query_seq, ref_seq, 10, 1, parasail.blosum62)
        #print(result.score1)
        #print(result.cigar)
        #print(result.ref_begin1)
        #print(result.ref_end1)
        #print(result.read_begin1)
        #print(result.read_end1)
        return result 
    
    def pairwise_sw_trace_align(self, ref_seq, query_seq):
        result = parasail.sw_trace(query_seq, ref_seq, 10, 1, parasail.blosum62)
        #print(result.cigar)
        #print(result.cigar.seq)
        print(result.cigar.decode)
        return result 
    
    def parse_cigar(self,result):
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
    
    def cigar_length(self,result):
        pairs = self.parse_cigar(result)
        total = 0
        for c,i in pairs:
            if c in ["M","=","X"]:
                total += i
            elif c in ["I","D","N","S","H","P"]:
                break
        return total
        
    def choose_closest_reference(self):
        closest_accession = None
        return closest_accession
    
    def identify_orf_coordinates(self, orf_coordinates, record_id=0):
        ref_sequence = self.get_reference_sequence(orf_coordinates, amino_acid=False)
        #print(len(ref_sequence))
        query_sequence = self.get_query_sequence(record_id, amino_acid=False)
        #print(len(query_sequence))
        result = self.pairwise_ssw_align(ref_sequence, query_sequence)
        return result.read_begin1, result.read_end1+1
    
    def discover_edits(self, orf_coordinates, found_coordinates, record_id=0):
        ref_sequence = self.get_reference_sequence(orf_coordinates)
        query_sequence = self.get_query_sequence(record_id, coordinates=found_coordinates)
        #print(ref_sequence[:30])
        #print(query_sequence[:30])
        #print(ref_sequence[-30:])
        #print(query_sequence[-30:])
        
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        cigar_length = self.cigar_length(result)
        while cigar_length < len(ref_sequence):
            print("cigar shorter than ref")
            found_coordinates, cigar_length, updated = self.frame_shift(orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_length, "", "N")
            if not updated:
                found_coordinates, cigar_length, updated = self.frame_shift(orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_length, "", "NN")
            if not updated:
                found_coordinates, cigar_length, updated = self.frame_shift(orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_length, "N", "")
            if not updated:
                break
        print(self.edits)
        return result
    
    def frame_shift(self, orf_coordinates, found_coordinates, record_id, ref_sequence, cigar_length, shift_from, shift_to):
        print("frame_shift from '", shift_from, "' to '", shift_to, "'")
        e = Edit(record_id, found_coordinates[0]+3*(cigar_length), shift_from, shift_to, self.closest_accession, orf_coordinates[0]+3*(cigar_length))
        e.apply_edit(self.consensus_sequence[record_id])
        coordinate_difference = len(shift_to) - len(shift_from)
        updated_found_coordinates = [found_coordinates[0], found_coordinates[1]+coordinate_difference]
        query_sequence = self.get_query_sequence(record_id, coordinates=updated_found_coordinates)
        result = self.pairwise_sw_trace_align(ref_sequence, query_sequence)
        new_cigar_length = self.cigar_length(result)
        updated = False
        
        #print(ref_sequence[-30:])
        #print(query_sequence[-30:])
        #print(cigar_length, len(ref_sequence))
        
        if new_cigar_length > cigar_length:
            updated = True
            self.edits.add_edit(e)
            return updated_found_coordinates, new_cigar_length, updated
        else:
            e.remove_edit(self.consensus_sequence[record_id])
            return found_coordinates, cigar_length, updated

    def run(self, reference_info_filepath, consensus_sequence_filepath, edit_filepath=""):
        self.load_input_files(consensus_sequence_filepath, reference_info_filepath, edit_filepath)
        print(self.reference_info["references"][self.closest_accession]["orf"])
        for key, value in self.reference_info["references"][self.closest_accession]["orf"].items():
            print(key,value)
            coordinates = (value["start"], value["end"])
            query_start, query_end = self.identify_orf_coordinates(orf_coordinates=coordinates)
            print(query_start, query_end)
            result = self.discover_edits(coordinates, (query_start, query_end))
