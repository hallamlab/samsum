__author__ = 'Connor Morgan-Lang'


import logging
import samsum
import numpy
import sys
import samsum.args as ss_args
import samsum.classy as ss_class
import samsum.logger as ss_log
import samsum.file_parsers as ss_fp
import samsum.utilities as ss_utils


def info(sys_args):
    """
    Function for writing version information about samsum and python dependencies.
    Other related info (citation, executable versions, etc.) should also be written through this sub-command.
    Create a SAMSumBase object for the `info` sub-command

    :param sys_args: List of arguments parsed from the command-line.
    :return: None
    """
    parser = ss_args.SAMSumArgumentParser(description="Return package and executable information.")
    args = parser.parse_args(sys_args)
    ss_log.prep_logging()
    info_ss = ss_class.SAMSumBase("info")

    logging.info("samsum version " + samsum.version + ".\n")

    # Write the version of all python deps
    py_deps = {"numpy": numpy.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # TODO: Write path to installation directory

    # Write the version of executable deps
    # info_ss.furnish_with_arguments(args)
    # logging.info(ss_utils.executable_dependency_versions(info_ss.executables))

    if args.verbose:
        pass
        # logging.info(summary_str)

    return


def stats(sys_args):
    """

    :param sys_args: List of arguments parsed from the command-line.
    :return: None
    """
    parser = ss_args.SAMSumArgumentParser(description="Calculate read coverage stats over reference sequences.")
    parser.add_stats_args()
    args = parser.parse_args(sys_args)
    # TODO: Create the log file based on the output table's path
    ss_log.prep_logging("samsum_log.txt", args.verbose)
    stats_ss = ss_class.SAMSumBase("stats")
    stats_ss.aln_file = args.am_file
    stats_ss.seq_file = args.fasta_file

    # Parse the FASTA file, calculating the length of each reference sequence and return this as a dictionary
    refseq_lengths = ss_fp.fasta_seq_lengths_ext(stats_ss.seq_file)
    references = ss_utils.load_references(refseq_lengths)
    refseq_lengths.clear()

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(stats_ss.aln_file, args.multireads)

    # TODO: Calculate the RPKM, FPKM and TPM for each reference sequence with reads mapped to it
    alignments = ss_utils.load_alignments(mapped_dict)
    mapped_dict.clear()
    ss_utils.load_reference_coverage(references, alignments)
    for ref_seq in references:  # type: ss_class.RefSequence
        print(ref_seq.get_info())
        sys.exit()

    # TODO: Calculate the percent sequence coverage for each reference sequence

    # TODO: Write the summary table with each of the above metrics as well as variance for each

    return
