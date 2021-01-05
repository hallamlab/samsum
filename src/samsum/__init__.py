"""
A light-weight python package for summarizing DNA sequence coverage from SAM files
"""

import _sam_module
from samsum import _version as samsum_version

name = "samsum"
__version__ = samsum_version.__version__
__all__ = ["_sam_module"]
