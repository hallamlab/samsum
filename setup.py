from glob import glob
from os.path import basename
from os.path import splitext
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

fasta_module = Extension("samsum._fasta_module",
                         sources=["src/extensions/fastamodule.cpp", "src/extensions/fastareader.cpp", "src/extensions/utilities.cpp"],
                         depends=["fastareader.h", "utilities.h", "types.h"],
                         include_dirs=["src/include/"],
                         language="c++")
sam_module = Extension("samsum._sam_module",
                       sources=["src/extensions/sammodule.cpp",
                                "src/extensions/helper.cpp", "src/extensions/sambamparser.cpp",
                                "src/extensions/utilities.cpp"],
                       depends=["helper.h", "sambamparser.h", "types.h", "utilities.h"],
                       include_dirs=["src/include/"],
                       language="c++")


SETUP_METADATA = \
    {
        "name": "samsum",
        "version": "0.0.8",
        "description": "A light-weight python package for summarizing sequence coverage from SAM and BAM files",
        "long_description": LONG_DESCRIPTION,
        "long_description_content_type": "text/markdown",
        "author": "Connor Morgan-Lang, Ryan McLaughlin",
        "author_email": "c.morganlang@gmail.com",
        "url": "https://github.com/hallamlab/samsum",
        "license": "GPL-3.0",
        "packages": find_packages('src', exclude=["tests"]),
        "include_package_data": True,
        "package_dir": {'samsum': 'src/samsum'},  # Necessary for proper importing
        "package_data": {'tests': ["tests/test-data/*.sam"]},
        "py_modules": [splitext(basename(path))[0] for path in glob('src/*.py')],
        "entry_points": {'console_scripts': ['samsum = samsum.__main__:main']},
        "classifiers": CLASSIFIERS,
        "ext_modules": [fasta_module, sam_module],
        "install_requires": ["numpy", "pytest"]
    }

setup(**SETUP_METADATA)
