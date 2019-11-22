__author__ = 'Connor Morgan-Lang'


class RefSequence:
    def __init__(self):
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
    def __init__(self, query_name: str, ref_name: str) -> None:
        self.query = query_name
        self.ref = ref_name
        self.start = 0
        self.end = 0
        self.percent_id = 0.0
        return

    def load_sam(self, aln_str: str) -> None:
        """
        Function for loading a string parsed from SAM alignments and is returned from _sam_module.get_mapped_reads
        :param aln_str: A string with various fields from a SAM file controlled by _sam_module.get_mapped_reads
        :return: None
        """
        fields = aln_str.split("\t")
        self.ref = fields[0]
        self.start = fields[1]
        return
