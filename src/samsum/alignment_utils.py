__author__ = 'Connor Morgan-Lang'

import logging
import sys
from samsum import classy


def load_references(refseq_lengths: dict) -> dict:
    """

    :param refseq_lengths:
    :return:
    """
    logging.debug("Loading the reference sequences into objects... ")
    references = {}
    for seq_name in refseq_lengths:  # type: str
        if seq_name in references:
            logging.error("Duplicate reference sequence names encountered: %s\n" % seq_name)
            sys.exit(3)
        ref_seq = classy.RefSequence(seq_name, refseq_lengths[seq_name])
        references[seq_name] = ref_seq
    logging.debug("done.\n")
    return references


def load_alignments(mapped_dict: dict) -> list:
    """

    :param mapped_dict:
    :return:
    """
    alignments = []
    for read_name in mapped_dict:
        for aln_dat in mapped_dict[read_name]:
            map_stats = aln_dat.split("\t")
            ref_seq = classy.AlignmentDat(read_name)
            ref_seq.load_sam(map_stats)
            alignments.append(ref_seq)
    return alignments


def load_reference_coverage(references: dict, alignments: list) -> None:
    """

    :param references:
    :param alignments:
    :return:
    """
    for aln in alignments:  # type: classy.AlignmentDat
        if aln.ref not in references:
            logging.error("Reference sequence from SAM file not found in FASTA: %s\n" % aln.ref)
            sys.exit(3)
        ref_seq = references[aln.ref]  # type: classy.RefSequence
        ref_seq.reads_mapped += 1
        ref_seq.weight_total += aln.weight
        if aln.start < ref_seq.leftmost:
            ref_seq.leftmost = aln.start
        if aln.end > ref_seq.rightmost:
            ref_seq.rightmost = aln.end
    return
