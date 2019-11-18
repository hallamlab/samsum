from setuptools import Extension
from setuptools import setup, find_packages

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

SETUP_METADATA = \
    {
        "name": "samsum",
        "version": "0.0.1",
        "description": "A light-weight python package for summarizing sequence coverage from SAM and BAM files",
        "long_description": LONG_DESCRIPTION,
        "long_description_content_type": "text/markdown",
        "author": "Connor Morgan-Lang, Kishori Konwar, Ryan McLaughlin",
        "author_email": "c.morganlang@gmail.com",
        "url": "https://github.com/hallamlab/samsum",
        "license": "GPL-3.0",
        "packages": find_packages(),
        "include_package_data": True,
        "entry_points": {'console_scripts': ['samsum = samsum.__main__:main']},
        "classifiers": CLASSIFIERS,
        "ext_modules": [Extension("_sam_parser",
                                  sources=["extensions/helper.c++", "extensions/matchoutputparser.c++",
                                           "extensions/rpkm.c++", "extensions/utilities.c++"],
                                  depends=["include/helper.h", "include/matchoutputparser.h",
                                           "include/rpkm.h", "include/types.h", "include/utilities.h"],
                                  include_dirs=["./samsum/include"],
                                  language="c++"),
                        Extension("_fasta_reader",
                                  sources=["extensions/fastareader.c++"],
                                  depends=["include/fastareader.h"],
                                  language="c++",
                                  include_dirs=["./samsum/include"])
                        ],
        "install_requires": ["numpy"]
    }

setup(**SETUP_METADATA)
