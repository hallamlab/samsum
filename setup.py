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

fasta_module = Extension("_fasta_module",
                         sources=["extensions/fastamodule.cpp", "extensions/fastareader.cpp"],
                         depends=["include/fastareader.h", "include/utilities.h", "types.h"],
                         include_dirs=["include/"],
                         language="c++")
sam_module = Extension("_sam_module",
                       sources=["extensions/sammodule.cpp",
                                "extensions/helper.cpp", "extensions/matchoutputparser.cpp",
                                "extensions/utilities.cpp"],
                       depends=["include/helper.h", "include/matchoutputparser.h",
                                "include/types.h", "include/utilities.h"],
                       include_dirs=["include/"],
                       language="c++")


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
        "ext_modules": [fasta_module, sam_module],
        "install_requires": ["numpy"]
    }

setup(**SETUP_METADATA)
