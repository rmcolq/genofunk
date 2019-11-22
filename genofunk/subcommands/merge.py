import logging
import os
import sys

from genofunk import merge


def run(options):
    if options.verbose:
        log_level = logging.DEBUG
        msg = "Using debug logging"
    else:
        log_level = logging.INFO
        msg = "Using info logging"

    #log_file = f"merge.log"
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
        "Input parameters:\nDirectory: %s\nOutput file: %s\nFeatures: %s\nMin occurence: %d\nInteractive: %s"
        % (options.directory, options.output_file, options.features, options.min_occurence, str(options.interactive)))

    g = merge.Merge()
    g.run(options.directory, options.output_file, options.features, options.min_occurence, options.interactive)
