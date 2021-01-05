# samsum
A light-weight python package for summarizing sequence coverage from SAM and BAM files

![tests](https://github.com/hallamlab/samsum/workflows/tests/badge.svg)
![build](https://github.com/hallamlab/samsum/workflows/build/badge.svg)
[![PyPI version](https://badge.fury.io/py/samsum.svg)](https://badge.fury.io/py/samsum)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4928c9ac353b4bdb93e351c0715a9fa1)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/samsum&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/hallamlab/samsum/branch/master/graph/badge.svg?token=q6UhWcnlS5)](https://codecov.io/gh/hallamlab/samsum)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/samsum/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/samsum/badges/platforms.svg)](https://anaconda.org/bioconda/samsum)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/samsum/badges/version.svg)](https://anaconda.org/bioconda/samsum)

## Installation

Samsum is currently supported on Mac and Linux systems and has been tested primarily on Ubuntu operating systems
 (bionic and trusty distributions).
It is a python package on the Python Package Index (PyPI) and can be installed using `pip`:

```pip install samsum```

Samsum can also be installed using `conda` with the command:

```bash
conda install -c bioconda samsum
```

You can also install samsum from source by cloning the directory from its GitHub page or downloading a GitHub release.

```bash
git clone https://github.com/hallamlab/samsum.git
cd samsum
python3 setup.py sdist
pip install dist/samsum*tar.gz
```

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
