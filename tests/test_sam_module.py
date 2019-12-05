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
    assert (len(mapping_list) == 14)
    return


def test_load_sam():
    from samsum import classy
    test_aln_data = ["refseq_name", "1", "1.0"]
    test_aln_obj = classy.AlignmentDat("readname")
    test_aln_obj.load_sam(test_aln_data)
    assert test_aln_obj.start == 1
    assert test_aln_obj.weight == 1.0
    assert test_aln_obj.ref == "refseq_name"
    return
