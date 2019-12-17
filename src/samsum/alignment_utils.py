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


def load_alignments(mapped_dict: dict) -> (list, int, int):
    """

    :param mapped_dict:
    :return:
    """
    alignments = []
    num_unmapped = 0
    mapped_sum = 0

    logging.info("Instantiating alignment data... ")

    for read_name in mapped_dict:
        for aln_dat in mapped_dict[read_name]:
            map_stats = aln_dat.split("\t")
            query_seq = classy.AlignmentDat(read_name)
            query_seq.load_sam(map_stats)
            if query_seq.ref == "UNMAPPED":
                num_unmapped = query_seq.weight
                continue
            alignments.append(query_seq)
            mapped_sum += query_seq.weight
    logging.info("done.\n")

    print("Num unique", len(mapped_dict))
    print("Mapped sum", mapped_sum)
    return alignments, num_unmapped, mapped_sum


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
        ref_seq.alignments.append(aln)
    return


def calculate_normalization_metrics(genome_dict: dict, sampled_reads: int) -> None:
    """
    Calculates the normalized abundance values for each header's RefSeq instance in genome_dict
        1. Reads per kilobase (RPK) is calculated using the reference sequence's length and number of reads (provided
        by the user via CLI)
        2. Fragments per kilobase per million mappable reads (FPKM) is calculated from the number of fragments
        (this is different from reads by, in a paired-end library, forward and reverse pair makes up one fragment)
        normalized by the reference sequence length and the number of reads mapped.
        2. Transcripts per million (TPM) is calculated similarly to FPKM but the order of operations is different.

    :param genome_dict: A dictionary of RefSeq instances indexed by headers (sequence names)
    :param sampled_reads: The number of reads sequenced (not aligned). Required for normalizing by sequencing depth
    :return: None
    """
    rpk_sum = 0  # The total reads per kilobase (RPK) of all reference sequences
    for header in sorted(genome_dict.keys()):
        ref_seq = genome_dict[header]  # type: classy.RefSequence
        if ref_seq.weight_total == 0:
            continue
        ref_seq.calc_fpkm(sampled_reads)
        ref_seq.rpk = sampled_reads/(ref_seq.length/1E3)
        rpk_sum += ref_seq.rpk

    denominator = rpk_sum / 1E6
    for header in genome_dict.keys():
        ref_seq = genome_dict[header]  # type: classy.RefSequence
        ref_seq.calc_tpm(denominator)

    return


def calculate_coverage(genome_dict: dict) -> None:
    for seq_name in genome_dict:  # type: str
        ref_seq = genome_dict[seq_name]  # type: classy.RefSequence
        if ref_seq.reads_mapped == 0:
            continue
        ref_seq.calc_coverage()
    return
