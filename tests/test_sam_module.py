import pytest
import os
import glob

from samsum import testing_utils as utils


@pytest.fixture()
def sam_list():
    data_path = utils.get_test_data("test-data")
    filenames = sorted(glob.glob(os.path.join(data_path, "*.sam")))
    return filenames


@pytest.fixture()
def alignment_dat_example():
    from samsum import classy
    ex_instance = classy.AlignmentDat("query_read_name")
    return ex_instance


@pytest.fixture()
def reference_sequence_example():
    from samsum import classy as ss_class
    refseq_ex = ss_class.RefSequence("NODE_1", 200)
    sam_dat = {"q1": ["NODE_1", "1", "5S45M", "0", "1.0"],
               "q2": ["NODE_1", "1", "45M", "0", "1.0"],
               "q3": ["NODE_1", "1", "5S145M", "0", "1.0"],
               "q4": ["NODE_1", "1", "5S145M", "0", "1.0"]}
    refseq_ex.reads_mapped = len(sam_dat)
    for r in sam_dat:
        aln_dat = ss_class.AlignmentDat(r)
        aln_dat.load_sam(sam_dat[r])
        refseq_ex.alignments.append(aln_dat)
    return refseq_ex


def test_get_mapped_reads():
    from samsum import _sam_module
    test_sam = utils.get_test_data('pytest_1.sam')
    assert isinstance(test_sam, str)
    mapping_list = _sam_module.get_mapped_reads(test_sam, False, 0)
    assert (len(mapping_list) == 16)
    return


def test_load_sam(alignment_dat_example):
    test_aln_data = ["refseq_name", "1", "5S145M", "0", "1.0"]
    alignment_dat_example.load_sam(test_aln_data)
    assert alignment_dat_example.start == 1
    assert alignment_dat_example.weight == 1.0
    assert alignment_dat_example.ref == "refseq_name"
    assert alignment_dat_example.end == 145
    assert alignment_dat_example.decode_cigar() == 145
    return


def test_load_unmapped(alignment_dat_example):
    test_unmapped_data = ['UNMAPPED', '-530976544', '', '0', '415858.0']
    alignment_dat_example.load_sam(test_unmapped_data)
    assert alignment_dat_example.weight == 415858.0
    assert alignment_dat_example.ref == "UNMAPPED"
    return


def test_decode_cigar(alignment_dat_example):
    cigar_str_1 = "101S19M30S"
    alignment_dat_example.cigar = cigar_str_1
    assert alignment_dat_example.decode_cigar() == 19
    assert alignment_dat_example.read_length == 150
    cigar_str_2 = "101S19M30H"
    alignment_dat_example.cigar = cigar_str_2
    assert alignment_dat_example.decode_cigar() == 19
    assert alignment_dat_example.read_length == 120

    return


def test_ref_sequence_abundances():
    from samsum import commands
    from samsum import classy as ss_class
    test_sam = utils.get_test_data("samsum_test_2.sam")
    test_aln = utils.get_test_data("samsum_test_2.fasta")
    ref_seq_abunds = commands.ref_sequence_abundances(aln_file=test_sam, seq_file=test_aln,
                                                      min_aln=10, p_cov=0, map_qual=0)
    e10 = ref_seq_abunds['AB-755_P22_E10_NODE_6_length_36342_cov_2156.57_ID_21']  # type: ss_class.RefSequence
    assert e10.reads_mapped == 10
    assert e10.weight_total == 5.0
    return


def test_proportion_covered(reference_sequence_example):
    assert reference_sequence_example.reads_mapped == 4
    proportion = reference_sequence_example.proportion_covered()
    assert proportion == 0.72

# TODO: Write test to ensure the RefSequence.rightmost doesn't exceed its length from alignment_dat_example.load_sam()


def test_output_table():
    from samsum import classy
    # TODO: Assert the header and the output rows contain the same number of fields
    # TODO: Assert the header follows the expected format
    return

