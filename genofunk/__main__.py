import argparse

import genofunk


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="genofunk",
        usage="genofunk <subcommand> <options>",
        description="Funky functions for virus consensus sequences",
    )

    parser.add_argument("--version", action="version", version=genofunk.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    #_____________________________ annotator ______________________________#
    subparser_annotator = subparsers.add_parser(
        "annotator",
        usage="genofunk annotator -r <reference_info_json> -f <consensus_fasta> [-a ref_accession -e <edit_file>]",
        help="Identify edits for a consensus sequence",
    )

    subparser_annotator.add_argument(
        "-r",
        "--reference_file",
        dest="reference_file",
        action="store",
        type=str,
        help="Input file: a JSON formatted file containing references. For each reference,"
        "the dictionary should include nucleotide sequence and ORF coordinates",
    )
    subparser_annotator.add_argument(
        "-f",
        "--consensus_file",
        dest="consensus_file",
        action="store",
        type=str,
        help="Input file: a fasta formatted file containing the consensus sequence(s) for "
        "a sample",
    )
    subparser_annotator.add_argument(
        "-a",
        "--accession",
        dest="accession",
        action="store",
        type=str,
        default=None,
        help="Accession (as written in the reference_file) for the closest reference sequence",
    )
    subparser_annotator.add_argument(
        "-e",
        "--edit_file",
        dest="edit_file",
        action="store",
        type=str,
        default=None,
        help="CSV file containing edits already found for this sample"
    )
    subparser_annotator.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    subparser_annotator.set_defaults(func=genofunk.subcommands.annotator.run)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
