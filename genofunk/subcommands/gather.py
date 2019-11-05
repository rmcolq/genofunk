import logging
import os
import sys

from genofunk import gather


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
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%d/%m/%Y %I:%M:%S",
    )
    logging.info(msg)
    logging.info(
        "Input parameters:\nDirectory: %s" %options.directory
    )

    g = gather.Gather()
    g.run(options.directory)
