import logging
import os
import sys

from genofunk import annotate, editfile


def run(options):
    if options.verbose:
        log_level = logging.DEBUG
        msg = "Using debug logging"
    else:
        log_level = logging.INFO
        msg = "Using info logging"

    #log_file = f"annotate.log"
    #if os.path.exists(log_file):
    #    os.unlink(log_file)
    logging.basicConfig(
        #filename=log_file,
        stream=sys.stdout,
        level=log_level,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%d/%m/%Y %I:%M:%S",
    )
    logging.info(msg)
    logging.info(
        "Input parameters:\nReference JSON: %s\nConsensus fasta: %s\nEdit file: %s\nAccession: %s\nStop codons: %s\n"
        "Min sequence length: %s\nNo stops in middle: %s\nFind compensating frame shifts: %s",
        options.reference_file,
        options.consensus_file,
        options.edit_file,
        options.accession,
        options.stop_codons,
        options.min_seq_length,
        options.no_stops_in_middle,
        options.find_compensating_frame_shifts
    )

    a = annotate.Annotate(options.accession)
    stop_codons = ",".split(options.stop_codons)
    a.run(options.reference_file, options.consensus_file, options.edit_file, stop_codons=stop_codons, max_mismatch=3,
          include_compensatory=options.find_compensating_frame_shifts, min_seq_length=options.min_seq_length,
          allow_stop_codons_in_middle = not options.no_stops_in_middle)
