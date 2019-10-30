import logging
import os
import sys

from genofunk import annotator, editfile


def run(options):
    if options.verbose:
        log_level = logging.DEBUG
        msg = "Using debug logging"
    else:
        log_level = logging.INFO
        msg = "Using info logging"

    #log_file = f"annotator.log"
    #if os.path.exists(log_file):
    #    os.unlink(log_file)
    logging.basicConfig(
        #filename=log_file,
        stream=sys.stdout,
        level=log_level,
        format="%(asctime)s %(message)s",
        datefmt="%d/%m/%Y %I:%M:%S",
    )
    logging.info(msg)
    logging.info(
        "Input parameters:\nReference JSON: %s\nConsensus fasta: %s\nEdit file: %s\nAccession: %s",
        options.reference_file,
        options.consensus_file,
        options.edit_file,
        options.accession,
    )

    a = annotator.Annotator(options.accession)
    a.run(options.reference_file, options.consensus_file, options.edit_file)

