__author__ = 'Connor Morgan-Lang'


class RefSequence:
    def __init__(self, ref_seq: str):
        self.name = ref_seq
        self.length = 0
        self.leftmost = 0
        self.rightmost = 0
        self.reads_mapped = 0
        self.weight_total = 0.0
        self.rpkm = 0.0
        self.fpkm = 0.0
        self.tpm = 0.0
        return

    def get_info(self):
        summary_str = "Reference sequence '%s':\n\t" % self.name
        summary_str += "\n\t".join(["Length = " + str(self.length) + "bp",
                                    "Number of reads mapped = " + str(self.reads_mapped),
                                    "Covered from %d to %d" % (self.leftmost, self.rightmost),
                                    "RPKM = %f" % self.rpkm,
                                    "FPKM = %f" % self.fpkm,
                                    "TPM  = %f" % self.tpm]) + "\n"
        return summary_str

    def calculate_coverage(self):
        return


class SAMSumBase:
    """
    A base class for all samsum sub-commands. It requires shared properties
    """
    def __init__(self, subcmd_name) -> None:
        self.subcmd = subcmd_name
        self.aln_file = ""
        self.seq_file = ""
        self.output_sep = ','
        return

    def get_info(self) -> str:
        info_string = ""
        return info_string


class AlignmentDat:
    """
    A class that stores alignmment information
    """
    def __init__(self, query_name: str) -> None:
        self.query = query_name
        self.ref = ""
        self.start = 0
        self.end = 0
        self.percent_id = 0.0
        self.weight = 0.0
        return

    def load_sam(self, aln_fields: list) -> None:
        """
        Function for loading a string parsed from SAM alignments and is returned from _sam_module.get_mapped_reads
        :param aln_fields: A list with various fields from a SAM file controlled by _sam_module.get_mapped_reads
        :return: None
        """
        fields = aln_fields
        self.ref = fields[0]
        self.start = fields[1]
        self.weight = fields[2]
        return
