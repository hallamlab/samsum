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


def load_alignments(mapped_dict: dict, min_aln: int) -> (list, float, float):
    """
    Converts the alignment strings for each query sequence into AlignmentDat instances. Sums the weights for unmapped
    (including those that fell below the minimum aligned percentage) and mapped reads.

    :param mapped_dict: A dictionary containing read names as keys and lists of SAM alignment rows as values. The format
     of these alignment rows is controlled by _sam_module.get_mapped_reads and must be accepted by AlignmentDat.load_sam
    :param min_aln: The minimum proportion of a read that must be aligned to a reference sequence to be included.
     If its aligned percentage falls below this threshold that query's alignment is not appended to the *alignments*
     list and its weight attribute is added to num_unmapped
    :return: A list of AlignmentDat instances, total alignment weights for unmapped reads and mapped reads
    """
    alignments = []
    num_unmapped = 0.0
    mapped_total = 0.0
    logging.info("Instantiating alignment data... ")

    for read_name in mapped_dict:
        for aln_dat in mapped_dict[read_name]:
            query_seq = classy.AlignmentDat(read_name, aln_dat.split("\t"))
            if query_seq.ref == "UNMAPPED":
                num_unmapped += query_seq.weight
                continue
            # Skip over alignments that don't meet or exceed the minimum percentage of the read length
            # TODO: Migrate this to the C++ side
            if 100*(query_seq.end-query_seq.start)/query_seq.read_length < min_aln:
                num_unmapped += query_seq.weight
                continue
            alignments.append(query_seq)
            mapped_total += query_seq.weight

    logging.info("done.\n")

    return alignments, num_unmapped, mapped_total


def load_reference_coverage(refseq_dict: dict, mapped_dict: dict, min_aln: int) -> (float, float):
    """
    Converts the alignment strings for each query sequence into AlignmentDat instances. Sums the weights for unmapped
    (including those that fell below the minimum aligned percentage) and mapped reads.

    :param refseq_dict: A dictionary of RefSequence instances indexed by headers (sequence names)
    :param mapped_dict: A dictionary containing read names as keys and lists of SAM alignment rows as values. The format
     of these alignment rows is controlled by _sam_module.get_mapped_reads and must be accepted by AlignmentDat.load_sam
    :param min_aln: The minimum proportion of a read that must be aligned to a reference sequence to be included.
     If its aligned percentage falls below this threshold that query's alignment is not appended to the *alignments*
     list and its weight attribute is added to num_unmapped
    :return: Total alignment weights for unmapped reads and mapped reads
    """
    logging.info("Associating read alignments with their respective reference sequences... ")
    num_unmapped = 0.0
    mapped_total = 0.0

    for refseq_name, alignment_data in mapped_dict.items():
        try:
            ref_seq = refseq_dict[refseq_name]  # type: classy.RefSequence
        except KeyError:
            if refseq_name != "UNMAPPED":
                logging.error("Reference sequence from SAM file not found in FASTA: %s\n" % refseq_name)
                sys.exit(3)
            else:
                unmapped_dat = classy.AlignmentDat(refseq_name, alignment_data.pop().split("\t"))
                num_unmapped += unmapped_dat.weight
                continue

        while alignment_data:  # type: list
            query_seq = classy.AlignmentDat(refseq_name, alignment_data.pop().split("\t"))

            if 100 * (query_seq.end - query_seq.start) / query_seq.read_length < min_aln:
                num_unmapped += query_seq.weight
                continue

            if query_seq.start < ref_seq.leftmost:
                ref_seq.leftmost = query_seq.start
            if query_seq.end > ref_seq.rightmost:
                ref_seq.rightmost = query_seq.end
            ref_seq.alignments.append(query_seq)

            ref_seq.reads_mapped += 1
            ref_seq.weight_total += query_seq.weight
            mapped_total += query_seq.weight
        ref_seq.calc_coverage()
        ref_seq.covered = ref_seq.proportion_covered()
        ref_seq.alignments.clear()

    logging.info("done.\n")
    return num_unmapped, mapped_total


def calculate_normalization_metrics(genome_dict: dict) -> None:
    """
    Calculates the normalized abundance values for each header's RefSeq instance in genome_dict
        1. Reads per kilobase (RPK) is calculated using the reference sequence's length and number of reads (provided
        by the user via CLI)
        2. Fragments per kilobase per million mappable reads (FPKM) is calculated from the number of fragments
        (this is different from reads by, in a paired-end library, forward and reverse pair makes up one fragment)
        normalized by the reference sequence length and the number of reads mapped.
        2. Transcripts per million (TPM) is calculated similarly to FPKM but the order of operations is different.

    :param genome_dict: A dictionary of RefSeq instances indexed by headers (sequence names)
    :return: None
    """
    rpk_sum = 0  # The total reads per kilobase (RPK) of all reference sequences
    for header in sorted(genome_dict.keys()):
        ref_seq = genome_dict[header]  # type: classy.RefSequence
        ref_seq.calc_rpk(ref_seq.weight_total)
        rpk_sum += ref_seq.rpk

    for header in sorted(genome_dict.keys()):
        ref_seq = genome_dict[header]  # type: classy.RefSequence
        if ref_seq.weight_total == 0:
            continue
        ref_seq.calc_fpkm(ref_seq.weight_total)

    denominator = rpk_sum / 1E6
    for header in genome_dict.keys():
        ref_seq = genome_dict[header]  # type: classy.RefSequence
        ref_seq.calc_tpm(denominator)

    return


def proportion_filter(references: dict, p_aln: int) -> float:
    """
    Removes all read alignments from a RefSequence with too little coverage, controlled by p_aln.
    The RefSequence.weight_total is then added to the discarded_weight, to be added to num_unmapped

    :param references: A dictionary of RefSequence instances indexed by headers (sequence names)
    :param p_aln: The minimum percentage of the reference sequence required to be covered by reads
    :return: The summed weight of each of the reads that were removed from low-coverage RefSequences
    """
    discarded_weight = 0.0
    logging.info("Filtering out reference sequences with coverage below " + str(p_aln) + "%... ")
    for seq_name in sorted(references):
        ref_seq = references[seq_name]  # type: classy.RefSequence
        if 100*ref_seq.covered < p_aln:
            discarded_weight += ref_seq.weight_total
            ref_seq.clear_alignments()
    logging.info("done.\n")
    return discarded_weight


def overlapping_intervals(coord_set_one, coord_set_two):
    s_one, e_one = coord_set_one
    s_two, e_two = coord_set_two
    if s_one <= s_two <= e_one or s_one <= e_two <= e_one:
        return True
    else:
        return False
