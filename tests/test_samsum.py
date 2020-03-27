import pytest


def test_overlapping_intervals():
    from samsum import alignment_utils
    coords_one = (0, 100)
    coords_two = (50, 101)
    coords_three = (101, 200)
    assert alignment_utils.overlapping_intervals(coords_one, coords_two) is True
    assert alignment_utils.overlapping_intervals(coords_two, coords_three) is True
    assert alignment_utils.overlapping_intervals(coords_one, coords_three) is False


@pytest.fixture()
def ref_sequence_abundances():
    from samsum import commands
    from samsum import testing_utils as utils
    test_sam = utils.get_test_data("samsum_test_2.sam")
    test_aln = utils.get_test_data("samsum_test_2.fasta")
    ref_seq_abunds = commands.ref_sequence_abundances(aln_file=test_sam, seq_file=test_aln,
                                                      min_aln=10, p_cov=0, map_qual=0)
    return ref_seq_abunds


def test_output_table(ref_sequence_abundances):
    from samsum import file_parsers as ss_fp
    curr_table_header = ["QueryName", "RefSequence", "ProportionCovered", "Coverage", "Fragments", "FPKM", "TPM"]
    ss_fp.write_summary_table(ref_sequence_abundances, "./samsum_table.csv", "pytest", 0)
    with open("samsum_table.csv") as table_handler:
        header_fields = table_handler.readline().strip().split(',')
        data_lines = []
        line = table_handler.readline()
        while line:
            data_lines.append(line.strip())
            line = table_handler.readline()

    assert header_fields == curr_table_header
    for line in data_lines:
        assert len(line.split(',')) == len(curr_table_header)
    return


def test_info():
    from samsum import __main__
    assert __main__.main(["samsum", "info"]) == 0
    assert __main__.main(["samsum"]) == 1
    return
