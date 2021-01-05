import os
import glob

import setuptools

package_root = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(package_root, "src", "samsum", "_version.py")) as fp:
    k, v = fp.read().strip().split(" = ")
version = v.strip('"')

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
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

extension = setuptools.Extension("_sam_module",
                                 sources=["src/extensions/sammodule.cpp",
                                          "src/extensions/helper.cpp", "src/extensions/sambamparser.cpp",
                                          "src/extensions/utilities.cpp", "src/extensions/types.cpp"],
                                 depends=["helper.h", "sambamparser.h", "types.h", "utilities.h"],
                                 include_dirs=["src/include/"],
                                 language="c++",
                                 extra_compile_args=[
                                     "-Wno-unused-result",
                                     "-Wno-cpp",
                                     "-Wno-unused-function",
                                 ]
                                 )


SETUP_METADATA = \
    {
        "name": "samsum",
        "version": version,
        "description": "A light-weight python package for summarizing sequence coverage from SAM and BAM files",
        "long_description": LONG_DESCRIPTION,
        "long_description_content_type": "text/markdown",
        "author": "Connor Morgan-Lang, Matthew Tang, Ryan J. McLaughlin",
        "author_email": "c.morganlang@gmail.com",
        "url": "https://github.com/hallamlab/samsum",
        "license": "GPL-3.0",
        "packages": setuptools.find_packages('src', exclude=["tests"]),
        "include_package_data": True,
        "package_dir": {'samsum': 'src/samsum'},  # Necessary for proper importing
        "package_data": {'tests': ["tests/test-data/*.sam"]},
        "py_modules": [os.path.splitext(os.path.basename(path))[0] for path in glob.glob('src/*.py')],
        "entry_points": {'console_scripts': ['samsum = samsum.__main__:main']},
        "classifiers": CLASSIFIERS,
        "ext_modules": [extension],
        "install_requires": ["numpy", "pytest", "pyfastx"]
    }

setuptools.setup(**SETUP_METADATA)
