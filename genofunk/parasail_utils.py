import logging
import parasail

def pairwise_sw_align(ref_seq, query_seq, gap_open=10, gap_extension=1, matrix=parasail.blosum62):
    """
    Run parasail-python SSW function
    :param ref_seq:
    :param query_seq:
    :param gap_open: Default 10
    :param gap_extension: Default 1
    :param matrix: Default BLOSUM62
    :return: parasail result (includes score, cigar, ref_begin, ref_end, read_begin, read_end attributes)
    """
    result = parasail.sw_trace_striped_sat(query_seq, ref_seq, gap_open, gap_extension, matrix)
    return result

def pairwise_nw_trace_align(ref_seq, query_seq, gap_open=10, gap_extension=1, matrix=parasail.blosum62):
    """
    :param ref_seq:
    :param query_seq:
    :param gap_open: Default 10
    :param gap_extension: Default 1
    :param matrix: Default BLOSUM62
    :return: parasail result (includes a cigar which can be decoded)
    """
    result = parasail.nw_trace(query_seq, ref_seq, gap_open, gap_extension, matrix)
    logging.debug("Parasail result %s" % result.cigar.decode)
    return result

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

def parse_cigar_pairs(result):
    """
    Extract the cigar from the parasail sw_trace_align or nw_trace_align result and turn it into cigar pairs
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

def get_alignment_start_end(result):
    # M   alignment match (can be a sequence match or mismatch)
    # I   insertion to the reference
    # D   deletion from the reference
    # N   skipped region from the reference
    # S   soft clipping (clipped sequences present in SEQ)
    # H   hard clipping (clipped sequences NOT present in SEQ)
    # P   padding (silent deletion from padded reference)
    # =   sequence match
    # X   sequence mismatch
    pairs = parse_cigar_pairs(result)
    ref_start = result.cigar.beg_ref
    ref_end = ref_start
    read_start = result.cigar.beg_query
    read_end = read_start
    found_alignment = False
    for c, i in pairs:
        if c in ["=", "X", "M"]:
            ref_end += i
            read_end += i
            found_alignment = True
        elif c in ["I"]:
            read_end += i
            if not found_alignment:
                read_start += i
        elif c in ["D", "N"]:
            ref_end += i
            if not found_alignment:
                ref_start += i
    return ref_start, ref_end - 1, read_start, read_end - 1

def get_cigar_length(pairs, max_mismatch, n_runs=[], min_match=3):
    """
    Find the length of aligned sequence until the first insertion/deletion/padding
    :param pairs:
    :return: number
    """
    #pairs = parse_cigar_pairs(result)
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

def cigar_has_no_indels(pairs):
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

def cigar_num_matches(pairs):
    """
    Are there non-match/mismatch symbols in cigar?
    :param pairs:
    :return: bool
    """
    total = 0
    for c, i in pairs:
        if c in ["="]:
            total += i
    return total

def cigar_score(pairs, match_score=1, mismatch_score=-1, gap_score=-1):
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

def cigar_has_min_matches(result, min_matches):
    pairs = parse_cigar_pairs(result)
    for c, i in pairs:
        if c in ["="] and i >= min_matches:
            return True
    return False


def case(c):
    #helper function
    if c in ["I", "D", "N", "S", "H", "P"]:
        return 1
    elif c in ["M", "X"]:
        return 2
    elif c == "=":
        return 3
    else:
        return 0

def is_extended_cigar_prefix(old_cigar_pairs, new_cigar_pairs):
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

def is_improved_cigar_prefix(old_cigar_pairs, new_cigar_pairs, consider_if_frameshift_added=True):
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
    if is_extended_cigar_prefix(old_cigar_pairs, new_cigar_pairs):
        logging.debug("Is extended cigar prefix, so an improvement")
        return True
    elif is_extended_cigar_prefix(new_cigar_pairs, old_cigar_pairs):
        logging.debug("Is cigar suffix, so NOT an improvement")
        return False

    old_cigar_num_matches = cigar_num_matches(old_cigar_pairs)
    new_cigar_num_matches = cigar_num_matches(new_cigar_pairs)

    # Second on existance of frameshift
    if consider_if_frameshift_added:
        old_cigar_has_no_indels = cigar_has_no_indels(old_cigar_pairs)
        new_cigar_has_no_indels = cigar_has_no_indels(new_cigar_pairs)
        logging.debug("Old cigar has no indels is %s, and number of matches is %i" %(old_cigar_has_no_indels, old_cigar_num_matches))
        logging.debug("New cigar has no indels is %s, and number of matches is %i" %(new_cigar_has_no_indels, new_cigar_num_matches))
        if old_cigar_has_no_indels and not new_cigar_has_no_indels \
                and old_cigar_num_matches >= new_cigar_num_matches:
            return False
        if new_cigar_has_no_indels and not old_cigar_has_no_indels \
                and new_cigar_num_matches >= old_cigar_num_matches:
            return True

    # Then on prefix comparison
    for i, (old_c, old_c_length) in enumerate(old_cigar_pairs):
        if i >= len(new_cigar_pairs):
            return False
        if new_cigar_pairs[i] == old_cigar_pairs[i]:
            continue

        new_c, new_c_length = new_cigar_pairs[i]

        if new_c == old_c and old_c_length - 1 <= new_c_length <= old_c_length:
            if new_cigar_num_matches > old_cigar_num_matches:
                return True
            elif old_cigar_num_matches > new_cigar_num_matches:
                return False

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
                                                                                 case(new_c) > case(old_c)))
        return case(new_c) > case(old_c)

    if len(new_cigar_pairs) > len(old_cigar_pairs):
        return True
    return False

def is_longer_cigar_prefix(old_cigar_pairs, new_cigar_pairs, max_mismatch):
    """
    Checks if new cigar prefix is longer than old cigar prefix. Only counts mismatch bases if there is a match after
    and breaks at first indel.
    :param old_cigar_pairs:
    :param new_cigar_pairs:
    :return: bool
    """
    old_cigar_length = get_cigar_length(old_cigar_pairs, max_mismatch)
    new_cigar_length = get_cigar_length(new_cigar_pairs, max_mismatch)
    logging.debug("Comparing %s and %s and found improvement to be %s" % (old_cigar_length, new_cigar_length,
                                                                              old_cigar_length < new_cigar_length))
    return old_cigar_length <= new_cigar_length