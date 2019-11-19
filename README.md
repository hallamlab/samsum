# samsum
A light-weight python package for summarizing sequence coverage from SAM and BAM files

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4928c9ac353b4bdb93e351c0715a9fa1)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/samsum&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/hallamlab/samsum/branch/master/graph/badge.svg?token=q6UhWcnlS5)](https://codecov.io/gh/hallamlab/samsum)

## Installation

Currently, it is not uploaded to the Python Package Index (PyPI) 
and needs to be built be being directed to the setup.py file like so:
```bash
python3 setup.py sdist
pip install dist/samsum*tar.gz
```


Once the package is available on PyPI it will be installed with the command
`pip install samsum`.

## Usage

`samsum stats` will read either a SAM or BAM file (this functionality will be implemented soon) and 
rapidly count the number of reads mapped to each reference sequence (e.g. contigs, scaffolds) 
while also keeping track of the reads that remain unmapped.
This all occurs within the C++ Python extension.
It will then read the reference FASTA file to gather the lengths of each reference sequence.
Combining the read counts and sequence lengths, it will then calculate:
- reads per kilobase per million mappable reads (RPKM)
- fragments per kilobase per million (FPKM)
- transcripts per milllion (TPM)

### Command-line options
By default, reads with multiple identical alignments (i.e. mapping quality is 0) are not included in these calculations.
This can be toggled off to include these alignments with the `--m` flag.
Another option is to drop counts for reference sequences if only a portion of a sequence is mapped to.
With the `-p` argument, you can control the minimum proportion a reference sequence needs to be covered 
for its read counts to be included in the output; all stats are otherwise set to 0.

An example command is:
```bash
samsum stats -c ref.fasta -a alignments.sam --m -p 0.5 -o output_dir/samsum_table.tsv
``` 

This will include all alignments, regardless of their mapping quality but only report alignments for reference sequences
that were covered across at least 50% of their length.

### API
 
Being a python package, samsum can also be readily imported into python code and used via its API.

.
.
.

## Outputs

The samsum_log.txt file is written to the output directory