__author__ = 'Connor Morgan-Lang'

import argparse


class SAMSumArgumentParser(argparse.ArgumentParser):
    """
    A base argparse ArgumentParser for samsum with functions to furnish with common arguments.
    This standardizes the interface for a unified aesthetic across all sub-commands
    """
    def __init__(self, **kwargs):
        """
        Instantiate the argparse argument-parser and create three broad argument groups:
            reqs - for the required parameters
            optopt - for the optional parameters
            miscellany - for the miscellaneous parameters that are module agnostic,
            for example verbose, help, num_threads
        :param kwargs:
        """
        super(SAMSumArgumentParser, self).__init__(add_help=False, **kwargs)
        self.reqs = self.add_argument_group("Required parameters")
        self.seqops = self.add_argument_group("Sequence operation arguments")
        self.optopt = self.add_argument_group("Optional options")
        self.miscellany = self.add_argument_group("Miscellaneous options")

        self.miscellany.add_argument("-v", "--verbose", action="store_true", default=False,
                                     help="Prints a more verbose runtime log")
        self.miscellany.add_argument("-h", "--help",
                                     action="help",
                                     help="Show this help message and exit")

    def parse_args(self, args=None, namespace=None):
        args = super(SAMSumArgumentParser, self).parse_args(args=args, namespace=namespace)

        return args

    def add_stats_args(self):
        self.reqs.add_argument("-f", "--ref_fasta",
                               required=True, dest="fasta_file",
                               help="Path to the reference file used to generate the SAM/BAM file.")
        self.reqs.add_argument("-a", "--alignments",
                               required=True, dest="am_file",
                               help="Path to a SAM/BAM file containing the read alignments to the reference FASTA.")
        self.seqops.add_argument("-l", "--aln_percent",
                                 required=False, dest="min_aln",
                                 default=10, type=int,
                                 help="The minimum percentage of a read's length that must be aligned to be included."
                                      " (DEFAULT = 10%%)")
        self.seqops.add_argument("-p", "--seq_coverage",
                                 required=False, dest="p_cov",
                                 default=50, type=int,
                                 help="The minimum percentage a reference sequence must be covered for its coverage"
                                      " stats to be included; they are set to zero otherwise. (DEFAULT = 50%%)")
        self.seqops.add_argument("-q", "--map_quality",
                                 required=False, dest="map_qual",
                                 default=0, type=int,
                                 help="The minimum mapping quality threshold for an alignment to pass. (DEFAULT = 0)")
        self.seqops.add_argument("--multireads",
                                 required=False,
                                 default=False, action="store_true",
                                 help="Flag indicating whether reads that mapped ambiguously to multiple positions"
                                      " (multireads) should be used in the counts.")
        self.optopt.add_argument("-o", "--output_table",
                                 required=False,
                                 default="./samsum_table.csv",
                                 help="Name of a file to write the alignment stats to."
                                      " (DEFAULT = ./samsum_table.csv)")
        self.optopt.add_argument("-s", "--sep",
                                 required=False,
                                 default=",", type=str,
                                 help="Field-separator character to be used when writing the output table."
                                      " (DEFAULT = ',')")
        return
