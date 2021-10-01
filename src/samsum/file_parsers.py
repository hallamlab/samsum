import os
import sys
import logging
import itertools

from pyfastx import Fasta
import pysam

import _sam_module
from samsum import classy as ss_class


LOGGER = logging.getLogger("samsum")


def open_pysam_alignment_file(aln_file: str, mode='r') -> pysam.AlignmentFile:
    file, ext = os.path.splitext(aln_file)
    if ext == '.bam':
        return pysam.AlignmentFile(aln_file, mode + 'b')
    elif ext == '.sam':
        return pysam.AlignmentFile(aln_file, mode)
    else:
        LOGGER.error("Unrecognized file extension in file '{}'\n".format(aln_file))
        sys.exit(11)


def pysam_parser(aln_file: str, multireads=False, aln_percent=0, min_mq=0) -> dict:
    sam_file = open_pysam_alignment_file(aln_file, 'r')
    unaligned_ref = "UNMAPPED"
    refseq_segments = {unaligned_ref: []}
    for segment in sam_file.fetch(until_eof=True):  # type: pysam.AlignedSegment
        if segment.reference_name:
            try:
                refseq_segments[segment.reference_name].append(segment)
            except KeyError:
                refseq_segments[segment.reference_name] = [segment]
        else:
            refseq_segments[unaligned_ref].append(segment)
    return refseq_segments


def sam_parser_ext(sam_file: str, multireads=False, aln_percent=0, min_mq=0) -> dict:
    """
    Wrapper function for using the _sam_parser extension to rapidly parse SAM files.

    :param sam_file: Path to the SAM file to be parsed
    :param multireads: Boolean flag indicating whether reads that have multiple ambiguous mapping positions are used
    :param aln_percent: The minimum percentage of a read's length that must be aligned to be included.
    :param min_mq: The minimum mapping quality for a read to be included in the analysis (as mapped)
    :return: A dictionary containing reference sequence names as keys and a list of Match instances as values
    """
    if not os.path.isfile(sam_file):
        LOGGER.error("SAM file '%s' doesn't exist.\n" % sam_file)
        sys.exit(3)

    reads_mapped = dict()
    mapping_list = iter(_sam_module.get_mapped_reads(sam_file, multireads, aln_percent, min_mq, 'r'))
    if not mapping_list:
        LOGGER.error("No alignments were read from SAM file '%s'\n" % sam_file)
        sys.exit(5)

    mapping_list_grouped = itertools.groupby(sorted(mapping_list, key=lambda x: x.subject), lambda x: x.subject)

    LOGGER.info("Grouping alignment data by reference sequence... ")

    for key, group in mapping_list_grouped:
        reads_mapped[key] = list(group)

    LOGGER.info("done.\n")

    LOGGER.debug("%d of unique read names returned by _sam_module.\n" % len(reads_mapped))

    return reads_mapped


def fasta_seq_lengths(fasta_file: str, min_seq_length=0) -> (dict, str):
    """
    Function for calculating the lengths of all sequences in a FASTA file.

    :param fasta_file: Path to a FASTA file to be parsed
    :param min_seq_length: The minimum length for a reference sequence to be included
    :return: A dictionary of sequence lengths indexed by their respective sequence names
    """
    seq_lengths_map = {}
    try:
        py_fa = Fasta(fasta_file, build_index=False, full_name=True)
    except RuntimeError as error:
        return seq_lengths_map, str(error)

    for name, seq in py_fa:  # type: (str, str)
        if len(seq) > min_seq_length:
            seq_lengths_map[name] = len(seq)

    return seq_lengths_map, ""


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

    ot_handler = open(output_table, 'w')

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
