import pytest
import samsum
import glob

from samsum import testing_utils as utils
from samsum import _sam_module

# @pytest.fixture()


def test_get_mapped_reads():
    test_sam = utils.get_test_data('pytest_1.sam')
    assert isinstance(test_sam, str)
    mapping_list = _sam_module.get_mapped_reads(test_sam, False, 0)
    assert (len(mapping_list) == 8)
