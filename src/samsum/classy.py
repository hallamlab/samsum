__author__ = 'Connor Morgan-Lang'

import logging
import sys
from samsum import utilities as ss_utils


class RefSequence:
    def __init__(self, ref_seq: str, seq_length: int):
        self.name = ref_seq
        self.length = seq_length
        self.leftmost = seq_length
        self.rightmost = 0
        self.reads_mapped = 0
        self.weight_total = 0.0
        self.rpk = 0.0
        self.fpkm = 0.0
        self.tpm = 0.0
        self.alignments = []
        return

    def get_info(self):
        summary_str = "Reference sequence '%s':\n\t" % self.name
        summary_str += "\n\t".join(["Length = " + str(self.length) + "bp",
                                    "Number of reads mapped = " + str(self.reads_mapped),
                                    "Covered from %d to %d" % (self.leftmost, self.rightmost),
                                    "RPKM = %f" % float(self.rpk/1E6),
                                    "FPKM = %f" % self.fpkm,
                                    "TPM  = %f" % self.tpm]) + "\n"
        return summary_str

    def aggregate(self, ref_seq) -> None:
        self.length += ref_seq.length
        self.leftmost = min([self.leftmost, ref_seq.leftmost])
        self.rightmost = min([self.rightmost, ref_seq.rightmost])
        self.weight_total += ref_seq.weight_total
        self.fpkm += ref_seq.fpkm
        self.rpk += ref_seq.rpk
        self.tpm += ref_seq.tpm
        self.alignments += ref_seq.alignments

    def proportion_covered(self) -> float:
        """
        Calculate the proportion of the RefSequence that was covered by mapped reads
        :return: Float representing the proportion of the Reference Sequence that was covered
        """
        return (self.rightmost-self.leftmost)/self.length

    def calc_coverage(self) -> float:
        """
        Calculate a 'dumb' coverage value of simply the number of base-pairs mapped
        (by summing the number of bp contained in the reads) and dividing by the length of sequence (RefSequence.length)
        :return:
        """
        bases_mapped = 0
        for aln_dat in self.alignments:  # type: AlignmentDat
            bases_mapped += (aln_dat.end - aln_dat.start)
        return bases_mapped/self.length

    def calc_rpk(self, num_reads):
        self.rpk = num_reads / (self.length / 1E3)

    def calc_fpkm(self, num_reads):
        mmr = float(num_reads/1E6)
        if self.weight_total == 0:
            self.fpkm = 0
        else:
            self.fpkm = float((self.weight_total/self.length)/mmr)
        return

    def calc_tpm(self, denominator) -> None:
        """
        Divide the read counts by the length of each gene in kilobases.
        Count up all the RPK values in a sample and divide this number by 1,000,000.
        Divide the RPK values by the “per million” scaling factor.

        :param denominator: The per-million scaling factor
        :return: None
        """
        if self.weight_total == 0:
            return
        self.tpm = self.rpk/denominator
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
        self.num_reads = 0
        return

    def get_info(self) -> str:
        info_string = "Info for " + self.subcmd + ":\n\t"
        info_string += "\n\t".join(["Alignment file: '%s'" % self.aln_file,
                                    "Reference sequence file: '%s'" % self.seq_file])
        if self.num_reads:
            info_string += "\n\tNumber of reads: %d" % self.num_reads
        return info_string + "\n"

    def furnish_with_arguments(self, args) -> None:
        self.executables["bwa"] = ss_utils.which("bwa")
        return


class AlignmentDat:
    """
    A class that stores alignment information
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
        consume_ref = ["M", "D", "N", "=", "X"]
        while i < len(self.cigar):
            if self.cigar[i].isdigit():
                buffer += self.cigar[i]
            elif buffer and self.cigar[i] in consume_ref:
                acc += int(buffer)
                buffer = ""
            else:
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
        self.end = self.start + self.cigar_length() - 1  # Need to subtract since SAM alignments are 1-based
        self.weight = float(fields[4])
        if self.weight > 1:
            logging.debug("Weight for '%s' is greater than 1 (%s).\n" % (self.query, str(self.weight)))
        return

    def get_info(self) -> str:
        info_string = "Info for alignment data:\n\t"
        info_string += "\n\t".join(["Query name: '%s'" % self.query,
                                    "Reference name: '%s'" % self.ref,
                                    "Start-End: %d - %d" % (self.start, self.end),
                                    "Weight: %d" % self.weight])
        return info_string
