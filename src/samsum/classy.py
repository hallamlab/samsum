__author__ = 'Connor Morgan-Lang'

import logging
from samsum import utilities as ss_utils
from samsum import alignment_utils as ss_aln_utils


class RefSequence:
    def __init__(self, ref_seq: str, seq_length: int):
        self.name = ref_seq
        self.length = seq_length
        self.leftmost = seq_length
        self.rightmost = 0
        self.reads_mapped = 0
        self.depth = 0.0
        self.covered = 0.0
        self.weight_total = 0.0
        self.rpk = 0.0
        self.fpkm = 0.0
        self.tpm = 0.0
        self.alignments = []
        self.tiles = []
        return

    def get_info(self):
        summary_str = "Reference sequence '%s':\n\t" % self.name
        summary_str += "\n\t".join(["Length = " + str(self.length) + "bp",
                                    "Number of reads mapped = " + str(self.reads_mapped) +
                                    " (%f fragments)" % self.weight_total,
                                    "Covered from %d to %d" % (self.leftmost, self.rightmost),
                                    "RPK = %f" % self.rpk,
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

    def merge_tiles(self) -> None:
        """
        Checks for Tile instances with overlapping ranges. Tile instances must have a 'start' and 'end' variable.

        :return: None
        """
        i = 0
        while i < len(self.tiles):
            tiles = sorted(self.tiles, key=lambda x: x.start)
            tile_i = tiles[i]  # type: Tile
            j = i + 1
            while j < len(tiles):
                tile_j = tiles[j]  # type: Tile
                if ss_aln_utils.overlapping_intervals((tile_i.start, tile_i.end), (tile_j.start, tile_j.end)):
                    # print("Merging:")
                    # print(tile_i.get_info(), "and", tile_j.get_info())
                    tile_i.merge(tile_j)
                    tiles.pop(j)
                else:
                    j += 1
            i += 1
        return

    def proportion_covered(self) -> float:
        """
        Calculate the proportion of the RefSequence that was covered by mapped reads.

        The algorithm works as follows:
            1. For each AlignmentDat instance in self.alignments:
                bin it into a continuously aligned regions (Tile)
            2. Merge the Tile instances from step one into the most contiguous possible
            3. Calculate the combined lengths of Tiles across the reference sequence and divide by its length

        :return: Float representing the proportion of the Reference Sequence that was covered
        """
        if self.reads_mapped == 0:
            return 0
        self.tiles.clear()
        for aln_dat in sorted(self.alignments, key=lambda x: x.start):  # type: AlignmentDat
            tile = Tile()
            tile.load_from_alignment_dat(aln_dat)
            i = 0
            while i < len(self.tiles):
                aln_coords = self.tiles[i]  # type: Tile
                if ss_aln_utils.overlapping_intervals((aln_coords.start, aln_coords.end), (aln_dat.start, aln_dat.end)):
                    aln_coords = self.tiles.pop(i)
                    tile.merge(aln_coords)
                    i = len(self.tiles)  # Increase i to the length of self.tiles to exit while loop
                i += 1
            self.tiles.append(tile)

        # Since coordinates are not compared as the tiles are merged the tiles can overlap necessitating a final merge
        self.merge_tiles()

        # Calculate the combined lengths of all tiles across the reference sequence
        total_tiled = 0
        for tile in self.tiles:
            total_tiled += (tile.end - tile.start)
        return total_tiled/self.length

    def calc_coverage(self) -> None:
        """
        Calculate a 'dumb' coverage value of simply the number of base-pairs mapped
        (by summing the number of bp contained in the reads) and dividing by the length of sequence (RefSequence.length)

        :return: None
        """
        bases_mapped = 0
        for aln_dat in self.alignments:  # type: AlignmentDat
            bases_mapped += (aln_dat.end - aln_dat.start)
        self.depth = bases_mapped/self.length
        return

    def calc_rpk(self, num_reads):
        self.rpk = num_reads / (self.length / 1E3)
        return

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

    def clear_alignments(self) -> None:
        self.reads_mapped = 0
        self.depth = 0.0
        self.weight_total = 0.0
        self.rpk = 0.0
        self.fpkm = 0.0
        self.tpm = 0.0
        self.alignments.clear()
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
        self.num_frags = 0
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


class Tile:
    def __init__(self):
        self.start = 0
        self.end = 0
        self.weight = 0.0

    def get_info(self):
        return "%d - %d: %f" % (self.start, self.end, self.weight)

    def load_from_alignment_dat(self, aln_dat):
        self.start = aln_dat.start
        self.end = aln_dat.end
        self.weight = aln_dat.weight

    def merge(self, another):
        self.start = min(self.start, another.start)
        self.end = max(self.end, another.end)
        self.weight += another.weight


class AlignmentDat(Tile):
    """
    A class that stores alignment information
    """
    def __init__(self, refseq_name: str, alignment_fields: list) -> None:
        super().__init__()
        self.ref = refseq_name
        self.query = ""
        self.cigar = ""
        self.read_length = 0
        self.percent_id = 0.0
        self.load_sam(alignment_fields)
        return

    def decode_cigar(self):
        self.read_length = 0
        aln_len = 0
        i = 0
        buffer = ""
        consume_ref = {"M", "D", "N", "=", "X"}  # type: set
        consume_query = {"M", "I", "S", "=", "X"}  # type: set
        while i < len(self.cigar):
            if self.cigar[i].isdigit():
                buffer += self.cigar[i]
            elif buffer:
                if self.cigar[i] in consume_ref:
                    aln_len += int(buffer)
                if self.cigar[i] in consume_query:
                    self.read_length += int(buffer)
                buffer = ""
            else:
                buffer = ""
            i += 1
        return aln_len

    def load_sam(self, aln_fields: list) -> None:
        """
        Function for loading a string parsed from SAM alignments and is returned from _sam_module.get_mapped_reads

        :param aln_fields: A list with various fields from a SAM file controlled by _sam_module.get_mapped_reads
        :return: None
        """
        fields = aln_fields
        self.query = fields[0]
        self.start = int(fields[1])
        self.cigar = fields[2]
        aln_len = self.decode_cigar()
        self.end = self.start + aln_len - 1  # Need to subtract since SAM alignments are 1-based
        self.weight = float(fields[4])
        if self.weight > 1:
            logging.debug("Weight for '%s' is greater than 1 (%s).\n" % (self.query, str(self.weight)))
        return

    def get_info(self) -> str:
        info_string = "Info for alignment data:\n\t"
        info_string += "\n\t".join(["Query name: '%s'" % self.query,
                                    "Reference name: '%s'" % self.ref,
                                    "Start-End: %d - %d" % (self.start, self.end),
                                    "Length: %d " % self.read_length,
                                    "Weight: %f" % self.weight])
        return info_string
