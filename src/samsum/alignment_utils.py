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


def load_alignments(mapped_dict: dict, min_aln: int) -> (list, int, int):
    """

    :param mapped_dict:
    :param min_aln:
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
                num_unmapped += query_seq.weight
                continue
            # Skip over alignments that don't meet or exceed the minimum percentage of the read length
            # TODO: Migrate this to the C++ side
            if 100*(query_seq.end-query_seq.start)/query_seq.read_length < min_aln:
                num_unmapped += query_seq.weight
                continue
            alignments.append(query_seq)
            mapped_sum += query_seq.weight
    logging.info("done.\n")

    return alignments, num_unmapped, mapped_sum


def load_reference_coverage(references: dict, alignments: list) -> None:
    """

    :param references:
    :param alignments:
    :return:
    """
    logging.info("Associating read alignments with their respective reference sequences... ")

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

    logging.info("done.\n")
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


def proportion_filter(references: dict, p_aln: int) -> float:
    """
    Removes all read alignments from a RefSequence with too little coverage, controlled by p_aln.
    The RefSequence.weight_total is then added to the discarded_weight, to be added to num_unmapped

    :param references: A dictionary of RefSeq instances indexed by headers (sequence names)
    :param p_aln: The minimum percentage of the reference sequence required to be covered by reads
    :return: The summed weight of each of the reads that were removed from low-coverage RefSequences
    """
    discarded_weight = 0.0
    logging.info("Filtering out reference sequences with coverage below " + str(p_aln) + "%... ")
    for seq_name in references:
        ref_seq = references[seq_name]  # type: classy.RefSequence
        if 100*ref_seq.proportion_covered() < p_aln:
            discarded_weight += ref_seq.weight_total
            ref_seq.clear_alignments()
    logging.info("done.\n")
    return discarded_weight
