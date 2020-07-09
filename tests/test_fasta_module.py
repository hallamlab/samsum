from samsum import testing_utils as utils


def test_ref_sequence_length() -> None:
    """
    Ensure the RefSequence.rightmost doesn't exceed its length from alignment_dat_example.load_sam()

    :return:
    """
    from samsum import file_parsers as ss_fp
    test_seq = utils.get_test_data("samsum_test_2.fasta")
    test_aln = utils.get_test_data("samsum_test_2.sam")
    ref_seq_lengths = ss_fp.fasta_seq_lengths_ext(fasta_file=test_seq)
    assert len(ref_seq_lengths) == 277

    # Parse the alignments and return the strings of reads mapped to each reference sequence
    mapped_dict = ss_fp.sam_parser_ext(test_aln, True)
    unmapped_aln = mapped_dict.pop("UNMAPPED").pop()  # type: _sam_module.Match
    assert unmapped_aln.weight == 4890.5

    # Convert the alignment strings returned by the sam_parser_ext into ss_class.AlignmentDat instances
    for ref_name, matches in mapped_dict.items():  # type: (str, list)
        for match in matches:
            assert match.end < ref_seq_lengths[match.subject]
    return


def test_fasta_reader() -> None:
    from samsum import file_parsers as ss_fp
    test_seq = utils.get_test_data("fasta_test.fa")
    ref_seq_lengths = ss_fp.fasta_seq_lengths_ext(fasta_file=test_seq)
    assert ref_seq_lengths == {"contig 1": 16, "contig 2": 25}
    return


if __name__ == "__main__":
    test_ref_sequence_length()
    test_fasta_reader()
