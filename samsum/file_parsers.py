__author__ = 'Connor Morgan-Lang'

import _sam_parser
import _fasta_reader
import os
import sys
import logging


def sam_parser_ext(sam_file: str, multireads=False) -> dict():
    """
    Wrapper function for using the _sam_parser extension to rapidly parse SAM files.

    :param sam_file: Path to the SAM file to be parsed
    :param multireads: Boolean flag indicating whether reads that have multiple ambiguous mapping positions are used
    :return: A dictionary mapping reference sequence names to the number of reads aligned to them
    """
    if not os.path.isfile(sam_file):
        logging.error("SAM file '%s' doesn't exist.\n" % sam_file)
        sys.exit(3)

    mapping_list = _sam_parser.refseq_read_counts(sam_file, multireads)
    if not mapping_list:
        logging.error("No alignments were read from SAM file '%s'\n" % sam_file)
        sys.exit(5)

    tmp_it = iter(mapping_list)
    reads_mapped = dict(zip(tmp_it, tmp_it))

    return reads_mapped


def fasta_seq_lengths_ext(fasta_file: str) -> dict():
    """
    Function for calculating the lengths of all sequences in a FASTA file.
    :param fasta_file: Path to a FASTA file to be parsed
    :return: A dictionary of sequence lengths indexed by their respective sequence names
    """
    if not os.path.isfile(fasta_file):
        logging.error("FASTA file '%s' doesn't exist.\n" % fasta_file)
        sys.exit(3)

    ext_seq_lengths = _fasta_reader.get_lengths(fasta_file)
    if not ext_seq_lengths:
        logging.error("No sequences were parsed from the FASTA file '%s'\n" % fasta_file)
        sys.exit(5)

    tmp_it = iter(ext_seq_lengths)
    seq_lengths_map = dict(zip(tmp_it, tmp_it))
    return seq_lengths_map
