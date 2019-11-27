import pytest
import samsum
import glob

from samsum import _sam_module

from . import testing_utils as utils


def test_get_mapped_reads():
    test_sam = utils.get_test_data('pytest_1.sam')
    test_sam = glob.glob(test_sam)
    mapping_list = _sam_module.get_mapped_reads(test_sam, False, 0.5)
    assert (len(mapping_list) == 16)
