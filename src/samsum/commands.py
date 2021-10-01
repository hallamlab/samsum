import os
import sys

from samsum import file_parsers as ss_fp
from samsum import alignment_utils as ss_aln_utils


def ref_sequence_abundances(aln_file: str, seq_file: str,
                            map_qual=0, p_cov=50, min_aln=10, multireads=False, logger=None,
                            **kwargs) -> dict:
    """
    An API function that will return a dictionary of RefSequence instances indexed by their sequence names/headers
    The RefSequence instances contain the populated variables:


    :param aln_file: Path to a SAM/BAM file containing the read alignments to the reference FASTA
    :param seq_file: Path to the reference FASTA file used to generate the SAM/BAM file
    :param map_qual: The minimum mapping quality threshold for an alignment to pass
    :param min_aln: The minimum percentage of a read's length that must be aligned to be included
    :param multireads: Flag indicating whether reads that mapped ambiguously to multiple positions (multireads)
    should be used in the counts
    :param p_cov: The minimum percentage a reference sequence must be covered for its coverage stats to be included;
    they are set to zero otherwise
    :param logger: A logging.Logger() instance named 'samsum'
    :param kwargs:
    :return: Dictionary of RefSequence instances indexed by their sequence names/headers
    """
    if not os.path.isfile(seq_file):
        logger.error("FASTA file '%s' doesn't exist.\n" % seq_file)
        sys.exit(3)

    logger.debug("Using Pyfastx to retrieve sequence lengths from FASTA... ")
    refseq_lengths, e = ss_fp.fasta_seq_lengths(seq_file)
    logger.debug("done.\n")

    if not refseq_lengths:
        logger.error("No sequences were parsed from the FASTA file '%s'\n" % seq_file)
        if e:
            logger.debug(e + "\n")
        sys.exit(5)
    else:
        logger.info(str(len(refseq_lengths)) + " sequences were read from " + seq_file + "\n")

    # Use these sequence lengths to load ReferenceSequence instances
    logger.debug("Loading the reference sequences into objects... ")
    references = ss_aln_utils.load_references(refseq_lengths)
    logger.debug("done.\n")
    refseq_lengths.clear()

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(aln_file, multireads, map_qual)

    logger.info("Associating read alignments with their respective reference sequences... ")
    num_unmapped, mapped_weight_sum = ss_aln_utils.load_reference_coverage(refseq_dict=references,
                                                                           mapped_dict=mapped_dict,
                                                                           min_aln=min_aln,
                                                                           log=logger)
    logger.info("done.\n")

    try:
        kwargs["mapped"] = mapped_weight_sum
        kwargs["num_frags"] = num_unmapped + mapped_weight_sum
    except KeyError:
        pass

    mapped_dict.clear()

    # Filter out alignments that with either short alignments or are from low-coverage reference sequences
    num_unmapped += ss_aln_utils.proportion_filter(references, p_cov)

    try:
        kwargs["unmapped"] = num_unmapped
    except KeyError:
        pass

    # Calculate the RPKM, FPKM and TPM for each reference sequence with reads mapped to it
    ss_aln_utils.calculate_normalization_metrics(references, num_unmapped)

    return references
