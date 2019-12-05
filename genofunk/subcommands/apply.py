import logging
import os
import sys

from genofunk import apply


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
    logging.info(msg)
    logging.info(
        "Input parameters:\nDirectory: %s\nEdit file: %s\nFeatures: %s\nOutput prefix: %s" %(options.directory,
                                                                                             options.edit_file,
                                                                                             options.features,
                                                                                             options.output_prefix)
    )

    g = apply.Apply()
    g.run(options.directory, options.edit_file, options.output_prefix, options.features)
