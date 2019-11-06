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

    #_____________________________ annotate ______________________________#
    subparser_annotate = subparsers.add_parser(
        "annotate",
        usage="genofunk annotate -r <reference_info_json> -f <consensus_fasta> -a ref_accession [-e <edit_file>]",
        help="Identify edits for a consensus sequence to remove frame shifts with respect to \n"
             "the reference. Assumes you know the closest reference for the provided sample \n"
             "and that this is the closest reference to all contigs in the consensus file",
    )

    subparser_annotate.add_argument(
        "-r",
        "--reference_file",
        dest="reference_file",
        action="store",
        type=str,
        help="Input file: a JSON formatted file containing references. For each reference,"
        "the dictionary should include nucleotide sequence and ORF coordinates",
    )
    subparser_annotate.add_argument(
        "-f",
        "--consensus_file",
        dest="consensus_file",
        action="store",
        type=str,
        help="Input file: a FASTA formatted file containing the consensus sequence(s) for "
        "a sample",
    )
    subparser_annotate.add_argument(
        "-a",
        "--accession",
        dest="accession",
        action="store",
        type=str,
        default=None,
        help="Accession (as written in the reference_file) for the closest reference sequence",
    )
    subparser_annotate.add_argument(
        "-e",
        "--edit_file",
        dest="edit_file",
        action="store",
        type=str,
        default=None,
        help="CSV file containing edits already found for this sample"
    )
    subparser_annotate.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    subparser_annotate.set_defaults(func=genofunk.subcommands.annotate.run)

    # _____________________________ merge ______________________________#
    subparser_merge = subparsers.add_parser(
        "merge",
        usage="genofunk merge -d <directory> ",
        help="Collates the edits for several files and flags frame shifts which occur multiple times and could be real",
    )

    subparser_merge.add_argument(
        "-d",
        "--directory",
        dest="directory",
        action="store",
        type=str,
        help="Input directory: containing pairs of consensus files and their edit lists (*.fasta, *.fasta.edits)",
    )

    subparser_merge.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    subparser_merge.set_defaults(func=genofunk.subcommands.merge.run)

    # _____________________________ apply ______________________________#
    subparser_apply = subparsers.add_parser(
        "apply",
        usage="genofunk apply -d <directory> -e <edit_file> -o <output_file>",
        help="Applies the edits to the consensus sequences and outputs an amino acid sequence for each",
    )

    subparser_apply.add_argument(
        "-d",
        "--directory",
        dest="directory",
        action="store",
        type=str,
        help="Input directory: containing pairs of consensus files and their edit lists (*.fasta, *.fasta.edits)",
    )
    subparser_apply.add_argument(
        "-e",
        "--edit_file",
        dest="edit_file",
        action="store",
        type=str,
        default=None,
        help="CSV file containing edits already found for this sample"
    )
    subparser_apply.add_argument(
        "-o",
        "--output",
        dest="output_file",
        action="store",
        type=str,
        help="Output file name",
    )

    subparser_apply.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    subparser_apply.set_defaults(func=genofunk.subcommands.apply.run)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
