# samsum
A light-weight python package for summarizing sequence coverage from SAM and BAM files

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4928c9ac353b4bdb93e351c0715a9fa1)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/samsum&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/hallamlab/samsum/branch/master/graph/badge.svg?token=q6UhWcnlS5)](https://codecov.io/gh/hallamlab/samsum)
[![Build Status](https://travis-ci.com/hallamlab/samsum.svg?branch=master)](https://travis-ci.com/hallamlab/samsum)

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

-   fragments per kilobase per million (FPKM)
-   transcripts per milllion (TPM)

### Command-line options
By default, reads with multiple identical alignments (i.e. mapping quality is 0) are not included in these calculations.
This can be toggled off to include these alignments with the `--m` flag.
Another option is to drop counts for reference sequences if only a portion of a sequence is mapped to.
With the `-p` argument, you can control the minimum proportion a reference sequence needs to be covered 
for its read counts to be included in the output; all stats are otherwise set to 0.

An example command is:
```bash
samsum stats -f ref.fasta -a alignments.sam --multireads -p 0.5 -o output_dir/samsum_table.tsv
``` 

This will include all alignments, regardless of their mapping quality but only report alignments for reference sequences
that were covered across at least 50% of their length.

### API
 
Being a python package, samsum can also be readily imported into python code and used via its API.

The function generally desired would be `ref_sequence_abundances`. Usage could be:
```python
from samsum import commands
sam="/home/user/reads_to_genome.sam"
fasta="/home/user/genome.fasta"
ref_seq_abunds = commands.ref_sequence_abundances(aln_file=sam, seq_file=fasta, min_aln=10, p_cov=0, map_qual=0)
```

The `ref_seq_abunds` object is a dictionary of `RefSequence` instances indexed by their header/sequence names.
`RefSequence` objects have several variables that are of interest:

-   `self.name` is the name of the (reference) sequence or header
-   `self.length` is the length (in base-pairs) of the sequence
-   `self.reads_mapped` is the number of reads that were mapped
-   `self.weight_total` is the number of fragments (float) that were mapped to the sequence
-   `self.fpkm` is Fragments Per Kilobase per Million mapped reads
-   `self.tpm` is Transcripts Per Million mapped reads

## Outputs

If `samsum stats` was executed, a "samsum_log.txt" file is written to the current working directory
 (i.e. where `samsum` was executed from). A comma-separated value (CSV) file with the fields 
 "QueryName", "RefSequence", "ProportionCovered", "Coverage", "Fragments", "FPKM" and "TPM" is written to a file
 path specified on the command-line, or by default "samsum_table.csv".
  A TSV file can be written instead if the `sep` argument was modified to 'tab'.
