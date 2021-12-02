from samsum import alignment as ss_aln_utils
from samsum import commands


class SAMSumBase(ss_aln_utils.Mapper):
    """
    A base class for all samsum sub-commands. It requires shared properties
    """

    def __init__(self, subcmd_name) -> None:
        super(SAMSumBase, self).__init__()
        self.subcmd = subcmd_name
        self.executables = {}
        self.ref_lengths = {}
        self.output_sep = ','

        # Alignment stats
        self.num_reads = 0
        self.num_frags = 0
        self.unmapped = 0
        self.mapped_weight = 0

        # Parsing options
        self.min_aln_percent = 10
        self.percent_coverage = 50

        return

    def get_info(self) -> None:
        info_string = "Info for " + self.subcmd + ":\n\t"
        info_string += "\n\t".join(["Alignment file: '%s'" % self.aln_file,
                                    "Reference sequence file: '%s'" % self.ref_fastx])
        if self.num_reads:
            info_string += "\n\tNumber of reads: %d" % self.num_reads

        self.log.debug(info_string + "\n")
        return

    def furnish_with_arguments(self, args) -> None:
        self.aln_file = args.am_file
        self.ref_fastx = args.fasta_file
        self.multireads = args.multireads
        self.min_mapping_q = args.map_qual

        return

    def calculate_abundances(self) -> dict:
        references = commands.ref_sequence_abundances(aln_file=self.aln_file, seq_file=self.ref_fastx,
                                                      map_qual=self.min_mapping_q, min_aln=self.min_aln_percent,
                                                      multireads=self.multireads, p_cov=self.percent_coverage,
                                                      logger=self.log,
                                                      num_frags=self.num_frags,
                                                      unmapped=self.unmapped, mapped=self.mapped_weight)

        self.get_info()

        return references
