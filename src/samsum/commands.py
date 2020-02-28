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


def ref_sequence_abundances(aln_file: str, seq_file: str, map_qual=0, p_cov=50, min_aln=10, multireads=False) -> dict:
    """
    An API function that will return a dictionary of RefSequence instances indexed by their sequence names/headers
    The RefSequence instances contain the populated variables:


    :param aln_file: Path to a SAM/BAM file containing the read alignments to the reference FASTA
    :param seq_file: Path to the reference FASTA file used to generate the SAM/BAM file
    :param map_qual: The minimum mapping quality threshold for an alignment to pass
    :param min_aln: The minimum percentage of a read's length that must be aligned to be included
    :param multireads: Flag indicating whether reads that mapped ambiguously to multiple positions (multireads) should be used in the counts
    :param p_cov: The minimum percentage a reference sequence must be covered for its coverage stats to be included; they are set to zero otherwise
    :return: Dictionary of RefSequence instances indexed by their sequence names/headers
    """
    refseq_lengths = ss_fp.fasta_seq_lengths_ext(seq_file)
    references = ss_aln_utils.load_references(refseq_lengths)
    refseq_lengths.clear()

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(aln_file, multireads, map_qual)

    ss_aln_utils.load_reference_coverage(refseq_dict=references, mapped_dict=mapped_dict, min_aln=min_aln)
    mapped_dict.clear()

    # Filter out alignments that with either short alignments or are from low-coverage reference sequences
    ss_aln_utils.proportion_filter(references, p_cov)

    # Calculate the RPKM, FPKM and TPM for each reference sequence with reads mapped to it
    ss_aln_utils.calculate_normalization_metrics(references)

    return references


def stats(sys_args):
    """
    A user-facing sub-command to write an abundance table from provided SAM and FASTA files.

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
    mapped_dict = ss_fp.sam_parser_ext(stats_ss.aln_file, args.multireads, args.map_qual)

    logging.debug(stats_ss.get_info())
    num_unmapped, mapped_weight_sum = ss_aln_utils.load_reference_coverage(refseq_dict=references,
                                                                           mapped_dict=mapped_dict,
                                                                           min_aln=args.min_aln)
    mapped_dict.clear()
    stats_ss.num_frags = num_unmapped + mapped_weight_sum

    # Filter out alignments that with either short alignments or are from low-coverage reference sequences
    num_unmapped += ss_aln_utils.proportion_filter(references, args.p_cov)

    # Calculate the RPKM, FPKM and TPM for each reference sequence with reads mapped to it
    ss_aln_utils.calculate_normalization_metrics(references)

    # Write the summary table with each of the above metrics as well as variance for each
    ss_fp.write_summary_table(references, args.output_table,
                              ss_utils.file_prefix(stats_ss.aln_file), num_unmapped, args.sep)

    # for seq_name in sorted(references):
    #     ref_seq = references[seq_name]  # type: ss_class.RefSequence
    #     if ref_seq.reads_mapped != 0:
    #         print(ref_seq.get_info())
    return
