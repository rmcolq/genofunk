import logging
import os
import sys

from genofunk import refparser


def run(options):
    if options.verbose:
        log_level = logging.DEBUG
        msg = "Using debug logging"
    else:
        log_level = logging.INFO
        msg = "Using info logging"

    #log_file = f"apply.log"
    #if os.path.exists(log_file):
    #    os.unlink(log_file)
    logging.basicConfig(
        #filename=log_file,
        stream=sys.stdout,
        level=log_level,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%d/%m/%Y %I:%M:%S",
    )
    if not options.output_file:
        options.output_file = options.genbank_file + ".json"
    logging.info(msg)
    logging.info(
        "Input parameters:\nGenbank file: %s\nOutput file: %s" %(options.genbank_file, options.output_file)
    )

    r = refparser.ReferenceParser()
    r.run(options.genbank_file, options.output_file)
