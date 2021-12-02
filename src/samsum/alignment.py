import os
import numpy as np
import sys
import logging
import psutil
from time import perf_counter
from datetime import timedelta

import mappy as mp
from tqdm import tqdm
from threading import Thread
from functools import partial

from samsum.multiprocessing import ThreadMap
from samsum import fastx_utils
from samsum.logger import CSVLogger

LOGGER = logging.getLogger("samsum")


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
                                    "FPKM = %f" % self.fpkm,
                                    "TPM  = %f" % self.tpm]) + "\n"
        return summary_str

    def aggregate(self, ref_seq) -> None:
        self.length += ref_seq.length
        self.leftmost = min([self.leftmost, ref_seq.leftmost])
        self.rightmost = min([self.rightmost, ref_seq.rightmost])
        self.weight_total += ref_seq.weight_total
        self.fpkm += ref_seq.fpkm
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
                if overlapping_intervals((tile_i.start, tile_i.end), (tile_j.start, tile_j.end)):
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
                if overlapping_intervals((aln_coords.start, aln_coords.end), (aln_dat.start, aln_dat.end)):
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
        return total_tiled / self.length

    def calc_coverage(self) -> None:
        """
        Calculate a 'dumb' coverage value of simply the number of base-pairs mapped
        (by summing the number of bp contained in the reads) and dividing by the length of sequence (RefSequence.length)

        :return: None
        """
        bases_mapped = 0
        for aln_dat in self.alignments:  # type: AlignmentDat
            bases_mapped += (aln_dat.end - aln_dat.start)
        self.depth = bases_mapped / self.length
        return

    def calc_fpkm(self, num_reads):
        mmr = float(num_reads / 1E6)
        if self.weight_total == 0:
            self.fpkm = 0
        else:
            self.fpkm = float((self.weight_total / self.length) / mmr)
        return

    def calc_tpm(self, denominator) -> None:
        """
        Divide the read counts by the length of each gene in kilobases.
        Count up all the FPK values in a sample and divide this number by 1,000,000.
        Divide the FPK values by the “per million” scaling factor.

        :param denominator: The per-million scaling factor
        :return: None
        """
        if self.weight_total == 0:
            return
        self.tpm = 1E6 * (self.fpkm / denominator)
        return

    def clear_alignments(self) -> None:
        self.reads_mapped = 0
        self.depth = 0.0
        self.weight_total = 0.0
        self.fpkm = 0.0
        self.tpm = 0.0
        self.alignments.clear()
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

    def decode_cigar(self) -> int:
        """
        Calculates and modifies the AlignmentDat read_length attribute from the cigar string attribute.
        The cigar string is parsed from the SAM file and represents the different states (e.g. insertion, deletion)
        and the length of these states across the aligned length.
        Also calculated in the alignment length itself and this length is returned as an integer.

        :return: An integer representing the alignment length
        """
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
        if self.weight > 1 and self.ref != "UNMAPPED":
            raise AssertionError("Weight for '%s' is greater than 1 (%s).\n" % (self.query, str(self.weight)))
        return

    def get_info(self) -> str:
        info_string = "Info for alignment data:\n\t"
        info_string += "\n\t".join(["Query name: '%s'" % self.query,
                                    "Reference name: '%s'" % self.ref,
                                    "Start-End: %d - %d" % (self.start, self.end),
                                    "Length: %d " % self.read_length,
                                    "Weight: %f" % self.weight])
        return info_string


class Mapper:
    """
    A class for performing sequence alignments
    """

    def __init__(self, num_threads=0):
        # Input parameters
        self.query_fastx = []  # Path to the query fastx file(s)
        self.ref_fastx = ""  # Path to the reference fasta file
        self.read_pairing = 2  # 2 == paired-end, 1 == single-end
        self.fq_fmt = 0  # 0 indicates separate files, 1 indicates interleaved paired-end reads
        self.threads = 0

        # Output parameters
        self.output_dir = ""  # Path for the indexed reference files and other temporary outputs
        self.log = LOGGER
        self.aln_f = ""  # Path to the alignment file (i.e. SAM or BAM)
        self.aln_format = "sam"

        # Alignment filtering options
        self.min_mapping_q = 0
        self.queries = 0
        self.matched_queries = 0
        self.multireads = False

        # Find maximum memory and set a RAM limit
        self.set_threads(num_threads)
        ram_mb = psutil.virtual_memory().total / (1024 * 1024)
        self.max_ram_perc = 20
        self.rlimit = ram_mb * (self.max_ram_perc / 100)
        self.log.info("System memory detected to be {} Mb."
                      "Setting RAM limit to {}% ({}).\n".format(ram_mb,
                                                                self.max_ram_perc,
                                                                self.rlimit))
        return

    def set_threads(self, num_threads: int) -> None:
        """Determine the number of threads to use based on system's cores."""
        _max_threads = 12
        # TODO: Find number of cores

        if self.threads > _max_threads:
            self.log.warning("Number of threads specified ({}) is greater than {}, which is not recommended."
                             "".format(self.threads, _max_threads))
        return

    def gen_pe_read_iter_fastx(self):
        return fastx_utils.gen_pe_read_iter_fastx(self.query_fastx, self.fq_fmt)

    def get_mappy_index(self):
        if self.multireads is False:
            max_hits = 1
        else:
            max_hits = 4
        return mp.Aligner(fn_idx_in=self.ref_fastx, best_n=max_hits, n_threads=self.threads)

    def idx_mappy_align(self) -> dict:
        """
        Single-threaded alignment with minimap2's Python binding, mappy.
        """
        hits = {}
        aln_idx = self.get_mappy_index()
        fwd_fq, rev_fq = self.gen_pe_read_iter_fastx()
        for fwd_read, rev_read in zip(fwd_fq, rev_fq):  # type: (tuple, tuple)
            self.queries += 1
            for hit in aln_idx.map(seq=fwd_read[1], seq2=rev_read[1]):  # type: mp.Alignment
                try:
                    hits[hit.ctg].append(hit)
                except KeyError:
                    hits[hit.ctg] = [hit]
            if self.queries % 1E4:
                if psutil.Process().memory_info().rss / (1024 * 1024) > self.rlimit:
                    LOGGER.error("RAM limit exceeded.\n")
                    sys.exit(1)

        self.report_completion_stats(hits)
        return hits

    def map_alignment_threads(self, aln_idx: mp.Aligner):
        """
        Align paired-end reads to a mappy Aligner object.
        """
        fwd_fq, rev_fq = self.gen_pe_read_iter_fastx()
        return ThreadMap(partial(MappyWorker, aln_idx),
                         zip(fwd_fq, rev_fq),
                         self.threads)

    def mappy_align_multi(self) -> dict:
        """
        Multi-processed implementation of mappy.
        :return:
        """
        aln_idx = self.get_mappy_index()

        results = self.map_alignment_threads(aln_idx)  # type: ThreadMap

        writer = AlignMapWriter(iterator=tqdm(results,
                                              desc="> calling",
                                              unit=" reads",
                                              leave=True),
                                aligner=aln_idx,
                                fd=self.aln_f)

        t0 = perf_counter()
        writer.start()
        writer.join()
        duration = perf_counter() - t0
        num_samples = sum(num_samples for read_id, num_samples in writer.log)

        self.report_completion_stats(writer.hits)
        return writer.hits

    def report_completion_stats(self, hits: dict):
        self.log.debug("Read-pairs queried: {}".format(self.queries))
        self.log.debug("Alignments: {}".format(sum([len(x) for x in hits.values()])))
        self.log.debug("Memory usage: {:.2f} Mb.".format(psutil.Process().memory_info().rss / (1024 * 1024)))
        # sys.stderr.write("> completed reads: %s\n" % len(writer.log))
        # sys.stderr.write("> duration: %s\n" % timedelta(seconds=np.round(duration)))
        # sys.stderr.write("> samples per second %.1E\n" % (num_samples / duration))
        # sys.stderr.write("> done\n")
        return


class AlignMapWriter(Thread):

    def __init__(self, iterator: iter, aligner: mp.Aligner, fd: str):
        super().__init__()
        if not fd:
            fd = sys.stdout
        self.fd = fd
        self.log = []
        self.aligner = aligner
        self.iterator = iterator
        self.write_headers()
        self.hits = {}
        return

    def write_headers(self):
        # if self.aligner:
        #     write_sam_header(self.aligner, fd=self.fd)
        pass

    def run(self):
        for read, hit in self.iterator:  # type: (str, mp.Alignment)
            if hit:
                try:
                    self.hits[hit.ctg].append(hit)
                except KeyError:
                    self.hits[hit.ctg] = [hit]
            # if len(seq):
            #     if self.aligner:
            #         print(read_id, seq)
            # write_sam(read_id, seq, qstring, mapping, fd=self.fd, unaligned=mapping is None)

            # self.log.append((read_id, samples))

            # else:
            #     LOGGER.warn("> skipping empty sequence %s", read_id)
        return


class MappyWorker(Thread):
    """
    Process that reads items from an input_queue, applies a func to them and puts them on an output_queue
    """

    def __init__(self, aligner, input_queue=None, output_queue=None):
        super().__init__()
        self.aligner = aligner
        self.input_queue = input_queue
        self.output_queue = output_queue

    def run(self):
        thrbuf = mp.ThreadBuffer()
        while True:
            item = self.input_queue.get()
            if item is StopIteration:
                self.output_queue.put(item)
                break
            fwd_mp, rev_mp = item
            mapping = next(self.aligner.map(seq=fwd_mp[1], seq2=rev_mp[1], buf=thrbuf), None)
            self.output_queue.put((fwd_mp[0], mapping))


def load_references(refseq_lengths: dict) -> dict:
    """

    :param refseq_lengths:
    :return:
    """
    references = {}
    for seq_name in refseq_lengths:  # type: str
        if seq_name in references:
            raise AssertionError("Duplicate reference sequence names encountered: %s\n" % seq_name)
        ref_seq = RefSequence(seq_name, refseq_lengths[seq_name])
        references[seq_name.split(' ')[0]] = ref_seq
    return references


def load_reference_coverage(refseq_dict: dict, mapped_dict: dict, min_aln: int, log) -> (float, float):
    """
    Converts the alignment strings for each query sequence into AlignmentDat instances. Sums the weights for unmapped
    (including those that fell below the minimum aligned percentage) and mapped reads.

    :param refseq_dict: A dictionary of RefSequence instances indexed by headers (sequence names)
    :param mapped_dict: A dictionary containing reference sequence names as keys and a list of Match instances as values
    :param min_aln: The minimum proportion of a read that must be aligned to a reference sequence to be included.
     If its aligned percentage falls below this threshold that query's alignment is not appended to the *alignments*
     list and its weight attribute is added to num_unmapped
    :param log: A logging.Logger instance to write messages to.
    :return: Total alignment weights for unmapped reads and mapped reads
    """
    num_unmapped = 0.0
    mapped_total = 0.0

    for refseq_name, alignment_data in mapped_dict.items():
        try:
            ref_seq = refseq_dict[refseq_name]  # type: RefSequence
        except KeyError:
            if refseq_name != "UNMAPPED":
                log.error("Reference sequence from SAM file not found in FASTA: %s\n" % refseq_name)
                sys.exit(3)
            else:
                unmapped_dat = alignment_data.pop()
                num_unmapped += unmapped_dat.weight
                continue

        while alignment_data:  # type: list
            query_seq = alignment_data.pop()

            if 100 * (query_seq.end - query_seq.start) / query_seq.read_length < min_aln:
                num_unmapped += query_seq.weight
                continue

            if query_seq.start < ref_seq.leftmost:
                ref_seq.leftmost = query_seq.start
            if query_seq.end > ref_seq.rightmost:
                ref_seq.rightmost = query_seq.end
            ref_seq.alignments.append(query_seq)

            ref_seq.reads_mapped += 1
            ref_seq.weight_total += query_seq.weight
            mapped_total += query_seq.weight
        ref_seq.calc_coverage()
        ref_seq.covered = ref_seq.proportion_covered()
        ref_seq.alignments.clear()

    return num_unmapped, mapped_total


def calculate_normalization_metrics(genome_dict: dict, unmapped_weight: float) -> None:
    """
    Calculates the normalized abundance values for each header's RefSeq instance in genome_dict
        1. Reads per kilobase (RPK) is calculated using the reference sequence's length and number of reads (provided
        by the user via CLI)
        2. Fragments per kilobase per million mappable reads (FPKM) is calculated from the number of fragments
        (this is different from reads by, in a paired-end library, forward and reverse pair makes up one fragment)
        normalized by the reference sequence length and the number of reads mapped.
        2. Transcripts per million (TPM) is calculated similarly to FPKM but the order of operations is different.

    :param genome_dict: A dictionary of RefSeq instances indexed by headers (sequence names)
    :param unmapped_weight: This represents the million-mappable reads for unaligned sequences. The 'weight' refers to
     this value being library-type agnostic; number of fragments (not reads!) for either a SE or PE library.
    :return: None
    """
    fpkm_sum = 0  # The total fragments per kilobase (FPK) of all reference sequences
    mill_frag_denom = unmapped_weight
    for header, ref_seq in sorted(genome_dict.items()):  # type: (str, RefSequence)
        mill_frag_denom += ref_seq.weight_total

    for header, ref_seq in sorted(genome_dict.items()):  # type: (str, RefSequence)
        if ref_seq.weight_total == 0:
            continue
        ref_seq.calc_fpkm(mill_frag_denom)
        fpkm_sum += ref_seq.fpkm

    for header in genome_dict.keys():
        ref_seq = genome_dict[header]  # type: RefSequence
        ref_seq.calc_tpm(fpkm_sum)

    return


def proportion_filter(references: dict, p_aln: int) -> float:
    """
    Removes all read alignments from a RefSequence with too little coverage, controlled by p_aln.
    The RefSequence.weight_total is then added to the discarded_weight, to be added to num_unmapped

    :param references: A dictionary of RefSequence instances indexed by headers (sequence names)
    :param p_aln: The minimum percentage of the reference sequence required to be covered by reads
    :return: The summed weight of each of the reads that were removed from low-coverage RefSequences
    """
    discarded_weight = 0.0
    logging.info("Filtering out reference sequences with coverage below " + str(p_aln) + "%... ")
    for seq_name in sorted(references):
        ref_seq = references[seq_name]  # type: RefSequence
        if 100 * ref_seq.covered < p_aln:
            discarded_weight += ref_seq.weight_total
            ref_seq.clear_alignments()
    logging.info("done.\n")
    return discarded_weight


def overlapping_intervals(coord_set_one, coord_set_two):
    s_one, e_one = coord_set_one
    s_two, e_two = coord_set_two
    if s_one <= s_two <= e_one or s_one <= e_two <= e_one:
        return True
    else:
        return False