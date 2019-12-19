__author__ = 'Connor Morgan-Lang'


import logging
import samsum
import numpy
import samsum.args as ss_args
import samsum.classy as ss_class
import samsum.logger as ss_log
import samsum.file_parsers as ss_fp
import samsum.utilities as ss_utils
import samsum.alignment_utils as ss_aln_utils


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
    info_ss.furnish_with_arguments(args)
    logging.info(ss_utils.executable_dependency_versions(info_ss.executables))

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
    references = ss_aln_utils.load_references(refseq_lengths)
    refseq_lengths.clear()

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(stats_ss.aln_file, args.multireads)

    # Convert the alignment strings returned by the sam_parser_ext into ss_class.AlignmentDat instances
    alignments, num_unmapped, mapped_weight_sum = ss_aln_utils.load_alignments(mapped_dict, args.min_aln)
    mapped_dict.clear()

    stats_ss.num_frags = num_unmapped + mapped_weight_sum
    logging.debug(stats_ss.get_info())
    ss_aln_utils.load_reference_coverage(references, alignments)
    alignments.clear()

    # Calculate the proportion sequence coverage for each reference sequence
    ss_aln_utils.calculate_coverage(references)

    # Filter out alignments that with either short alignments or are from low-coverage reference sequences
    num_unmapped += ss_aln_utils.proportion_filter(references, args.p_cov)

    # Calculate the RPKM, FPKM and TPM for each reference sequence with reads mapped to it
    ss_aln_utils.calculate_normalization_metrics(references, stats_ss.num_frags)

    # Write the summary table with each of the above metrics as well as variance for each
    ss_fp.write_summary_table(references, args.output_table,
                              ss_utils.file_prefix(stats_ss.aln_file), num_unmapped, args.sep)

    return
