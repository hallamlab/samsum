"""Modified from sourmash's sourmash/tests/sourmash_tst_utils.py"""

import os
from pkg_resources import Requirement, resource_filename, ResolutionError


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("samsum"), "tests/test-data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'tests/test-data', filename)
    return filepath


def get_project_root():
    return resource_filename(Requirement.parse("samsum"), "")
