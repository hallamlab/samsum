import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from samsum import alignment as aln_utils
        self.pe_query_fastx = [get_test_data("samsum_test_2.R1.fq"), get_test_data("samsum_test_2.R2.fq")]
        self.il_query_fastx = [get_test_data("SI072_40K.fastq")]
        # self.il_query_fastx = ["/media/connor/Rufus/Sakinaw/RawData/Joint_Genome_Institute/metagenomes/sak_2013_06_06_120m/12414.4.257161.TGACTGA-GTCAGTC.filter-METAGENOME.fastq.gz"]
        # self.il_query_fastx = ["/media/connor/Rufus/Sakinaw/RawData/Joint_Genome_Institute/metatranscriptomes/sak_2013_06_06_120m/12600.5.269981.ATGACGT-GACGTCA.filter-MTF.fastq.gz"]

        self.mapper = aln_utils.Mapper()
        self.mapper.ref_fastx = get_test_data("samsum_test_2.fasta")

        return

    def test_overlapping_intervals(self):
        from samsum import alignment
        coords_one = (0, 100)
        coords_two = (50, 101)
        coords_three = (101, 200)
        self.assertTrue(alignment.overlapping_intervals(coords_one, coords_two))
        self.assertTrue(alignment.overlapping_intervals(coords_two, coords_three))
        self.assertFalse(alignment.overlapping_intervals(coords_one, coords_three))
        return

    def test_mappy_align(self):
        self.mapper.query_fastx = self.il_query_fastx
        self.mapper.fq_fmt = 1
        self.mapper.mappy_align(num_threads=8)
        return

    def test_idx_mappy_align(self):
        self.mapper.query_fastx = self.il_query_fastx
        self.mapper.fq_fmt = 1
        ref_aln_counts = self.mapper.idx_mappy_align()
        self.assertEqual(277, len(ref_aln_counts))
        self.assertEqual(2000, sum(ref_aln_counts.values()))
        return

    def test_get_pe_alignments(self):
        from samsum import fastx_utils
        self.mapper.query_fastx = self.pe_query_fastx
        ref_fasta = fastx_utils.gen_fastx(self.mapper.ref_fastx)
        total_hits = 0
        for name, seq in ref_fasta:
            hits = self.mapper.get_pe_alignments(ref_seq=seq, ref_name=name)
            total_hits += hits[name]
        self.assertEqual(2000, total_hits)

        # Test interleaved fastq
        self.mapper.query_fastx = self.il_query_fastx
        self.mapper.fq_fmt = 1
        for name, seq in ref_fasta:
            hits = self.mapper.get_pe_alignments(ref_seq=seq, ref_name=name)
            total_hits += hits[name]
        self.assertEqual(2106, total_hits)
        return


if __name__ == '__main__':
    unittest.main()
