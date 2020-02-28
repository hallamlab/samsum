
import os
import sys
import logging
from samsum import _fasta_module, _sam_module
from samsum import classy as ss_class

__author__ = 'Connor Morgan-Lang'


def sam_parser_ext(sam_file: str, multireads=False, aln_percent=0, min_mq=0) -> dict:
    """
    Wrapper function for using the _sam_parser extension to rapidly parse SAM files.

    :param sam_file: Path to the SAM file to be parsed
    :param multireads: Boolean flag indicating whether reads that have multiple ambiguous mapping positions are used
    :param aln_percent: The minimum percentage of a read's length that must be aligned to be included.
    :param min_mq: The minimum mapping quality for a read to be included in the analysis (as mapped)
    :return: A dictionary mapping query sequence (read) names to a list of alignment data strings
    """
    if not os.path.isfile(sam_file):
        logging.error("SAM file '%s' doesn't exist.\n" % sam_file)
        sys.exit(3)

    reads_mapped = dict()
    mapping_list = iter(_sam_module.get_mapped_reads(sam_file, multireads, aln_percent, min_mq, 'r'))
    if not mapping_list:
        logging.error("No alignments were read from SAM file '%s'\n" % sam_file)
        sys.exit(5)

    logging.info("Zipping query names and alignment data... ")
    for ref, aln_dat in zip(mapping_list, mapping_list):
        try:
            reads_mapped[ref].append(aln_dat)
        except KeyError:
            reads_mapped[ref] = [aln_dat]
    logging.info("done.\n")

    logging.debug("%d of unique read names returned by _sam_module.\n" % len(reads_mapped))

    return reads_mapped


def fasta_seq_lengths_ext(fasta_file: str, min_seq_length=0) -> dict:
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


def write_summary_table(references: dict, output_table: str, samsum_exp: str, unmapped_reads: float, sep=",") -> None:
    """
    Writes the output file most people care about - the table summarizing abundance metrics for each reference sequence.
    Takes a dictionary of sequence names indexing their RefSequence instances and writes specific data for each.
    Current header is:
    [RefSequence.name, Query.name, ProportionCovered, Reads, RPKM, FPKM, TPM]
    Included in this table as the first row are the unmapped reads (UNMAPPED) with relevant information where possible

    :param references: A dictionary of RefSequence instances indexed by the reference sequence names (headers)
    :param samsum_exp: String representing the origin of the query reads, or alignment experiment name
    :param output_table: A string representing the path of the file to write to
    :param unmapped_reads: The number of reads that were not mapped to the reference sequences
    :param sep: Field separator to use. The default is a comma.
    :return: None
    """
    header = ["QueryName", "RefSequence", "ProportionCovered", "Coverage", "Fragments", "FPKM", "TPM"]
    buffer = sep.join(header) + "\n"
    # Add the unmapped reads data
    buffer += sep.join([samsum_exp, "UNMAPPED", "NA", "NA", str(unmapped_reads), "NA", "NA"]) + "\n"

    try:
        ot_handler = open(output_table, 'w')
    except IOError:
        logging.error("Unable to open output table '%s' for writing.\n" % output_table)
        sys.exit(3)

    for seq_name in sorted(references, key=lambda x: references[x].tpm, reverse=True):  # type: str
        ref_seq = references[seq_name]  # type: ss_class.RefSequence
        data_fields = [ref_seq.covered, ref_seq.depth, ref_seq.weight_total, ref_seq.fpkm, ref_seq.tpm]

        buffer += sep.join([samsum_exp, ref_seq.name] +
                           [str(round(x, 3)) for x in data_fields]) + "\n"
        if len(buffer) > 1E6:
            ot_handler.write(buffer)
            buffer = ""
    ot_handler.write(buffer)

    return
