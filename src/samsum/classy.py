__author__ = 'Connor Morgan-Lang'


from samsum import utilities as ss_utils


class RefSequence:
    def __init__(self, ref_seq: str, seq_length: int):
        self.name = ref_seq
        self.length = seq_length
        self.leftmost = seq_length
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
        self.executables = {}
        self.aln_file = ""
        self.seq_file = ""
        self.output_sep = ','
        return

    def get_info(self) -> str:
        info_string = ""
        return info_string

    def furnish_with_arguments(self, args) -> None:
        self.executables["bwa"] = ss_utils.which("bwa")
        return


class AlignmentDat:
    """
    A class that stores alignmment information
    """
    def __init__(self, query_name: str) -> None:
        self.query = query_name
        self.ref = ""
        self.cigar = ""
        self.start = 0
        self.end = 0
        self.percent_id = 0.0
        self.weight = 0.0
        return

    def cigar_length(self):
        acc = 0
        i = 0
        buffer = ""
        while i < len(self.cigar):
            if self.cigar[i].isdigit():
                buffer += self.cigar[i]
            elif buffer:
                acc += int(buffer)
                buffer = ""
            i += 1
        return acc

    def load_sam(self, aln_fields: list) -> None:
        """
        Function for loading a string parsed from SAM alignments and is returned from _sam_module.get_mapped_reads
        :param aln_fields: A list with various fields from a SAM file controlled by _sam_module.get_mapped_reads
        :return: None
        """
        fields = aln_fields
        self.ref = fields[0]
        self.start = int(fields[1])
        self.cigar = fields[2]
        self.end = self.start + self.cigar_length()
        self.weight = float(fields[3])
        return
