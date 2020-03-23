from samsum import testing_utils as utils


def test_ref_sequence_length():
    """
    Ensure the RefSequence.rightmost doesn't exceed its length from alignment_dat_example.load_sam()

    :return:
    """
    from samsum import file_parsers as ss_fp
    from samsum import alignment_utils as ss_aln_utils
    from samsum import classy
    test_seq = utils.get_test_data("samsum_test_2.fasta")
    test_aln = utils.get_test_data("samsum_test_2.sam")
    ref_seq_lengths = ss_fp.fasta_seq_lengths_ext(fasta_file=test_seq)
    assert len(ref_seq_lengths) == 277

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(test_aln, True)
    # Convert the alignment strings returned by the sam_parser_ext into ss_class.AlignmentDat instances
    alignments, num_unmapped, mapped_weight_sum = ss_aln_utils.load_alignments(mapped_dict, 10)
    for aln_dat in alignments:  # type: classy.AlignmentDat
        assert aln_dat.end < ref_seq_lengths[aln_dat.ref]
    return


def test_fasta_reader():
    from samsum import file_parsers as ss_fp
    test_seq = utils.get_test_data("fasta_test.fa")
    ref_seq_lengths = ss_fp.fasta_seq_lengths_ext(fasta_file=test_seq)
    assert ref_seq_lengths == {"contig_1": 16, "contig_2": 25}
    return
