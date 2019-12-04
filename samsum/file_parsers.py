
import os
import sys
import logging
from samsum import _fasta_module, _sam_module

__author__ = 'Connor Morgan-Lang'


def sam_parser_ext(sam_file: str, multireads=False, min_mq=0) -> dict():
    """
    Wrapper function for using the _sam_parser extension to rapidly parse SAM files.

    :param sam_file: Path to the SAM file to be parsed
    :param multireads: Boolean flag indicating whether reads that have multiple ambiguous mapping positions are used
    :param min_mq: The minimum mapping quality for a read to be included in the analysis (as mapped)
    :return: A dictionary mapping reference sequence names to the number of reads aligned to them
    """
    if not os.path.isfile(sam_file):
        logging.error("SAM file '%s' doesn't exist.\n" % sam_file)
        sys.exit(3)

    # TODO: Convert to a generator function so _sam_parser returns chunks in a list of 10,000 read strings
    # mapping_list = _sam_parser.refseq_read_counts(sam_file, multireads)
    mapping_list = _sam_module.get_mapped_reads(sam_file, multireads, min_mq)
    if not mapping_list:
        logging.error("No alignments were read from SAM file '%s'\n" % sam_file)
        sys.exit(5)

    tmp_it = iter(mapping_list)
    reads_mapped = dict(zip(tmp_it, tmp_it))

    return reads_mapped


def fasta_seq_lengths_ext(fasta_file: str, min_seq_length=0) -> dict():
    """
    Function for calculating the lengths of all sequences in a FASTA file.
    :param fasta_file: Path to a FASTA file to be parsed
    :param min_seq_length: The minimum length for a reference sequence to be included
    :return: A dictionary of sequence lengths indexed by their respective sequence names
    """
    if not os.path.isfile(fasta_file):
        logging.error("FASTA file '%s' doesn't exist.\n" % fasta_file)
        sys.exit(3)

    logging.debug("Using FASTA module to retrieve sequence lengths from FASTA... ")
    ext_seq_lengths = _fasta_module.get_lengths(fasta_file, min_seq_length)
    if not ext_seq_lengths:
        logging.error("No sequences were parsed from the FASTA file '%s'\n" % fasta_file)
        sys.exit(5)
    logging.debug("done.\n")

    logging.debug("Converting list of sequence lengths into dictionary... ")
    tmp_it = iter(ext_seq_lengths)
    seq_lengths_map = dict(zip(tmp_it, tmp_it))
    logging.debug("done.\n")

    logging.info(str(len(seq_lengths_map)) + " sequences were read from " + fasta_file + "\n")

    return seq_lengths_map
