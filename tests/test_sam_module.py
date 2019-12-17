import pytest
import os
import glob

from samsum import testing_utils as utils


@pytest.fixture()
def sam_list():
    data_path = utils.get_test_data("test-data")
    filenames = sorted(glob.glob(os.path.join(data_path, "*.sam")))
    return filenames


def test_get_mapped_reads():
    from samsum import _sam_module
    test_sam = utils.get_test_data('pytest_1.sam')
    assert isinstance(test_sam, str)
    mapping_list = _sam_module.get_mapped_reads(test_sam, False, 0)
    assert (len(mapping_list) == 16)
    return


@pytest.fixture()
def alignment_dat_example():
    from samsum import classy
    ex_instance = classy.AlignmentDat("query_read_name")
    return ex_instance


def test_load_sam(alignment_dat_example):
    test_aln_data = ["refseq_name", "1", "5S145M", "0", "1.0"]
    alignment_dat_example.load_sam(test_aln_data)
    assert alignment_dat_example.start == 1
    assert alignment_dat_example.weight == 1.0
    assert alignment_dat_example.ref == "refseq_name"
    return


def test_load_unmapped(alignment_dat_example):
    test_unmapped_data = ['UNMAPPED', '-530976544', '', '0', '415858.0']
    alignment_dat_example.load_sam(test_unmapped_data)
    assert alignment_dat_example.weight == 415858.0
    assert alignment_dat_example.ref == "UNMAPPED"
    return


def test_cigar_length(alignment_dat_example):
    cigar_str = "101S19M30S"
    alignment_dat_example.cigar = cigar_str
    assert alignment_dat_example.cigar_length() == 19
    return


# TODO: Write test to ensure the RefSequence.rightmost doesn't exceed its length from alignment_dat_example.load_sam()
