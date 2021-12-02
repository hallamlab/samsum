import os
import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from samsum import alignment as aln_utils
        self.pe_query_fastx = [get_test_data("samsum_test_2.R1.fq"), get_test_data("samsum_test_2.R2.fq")]
        self.il_query_fastx = [get_test_data("SI072_40K.fastq")]

        self.mapper = aln_utils.Mapper()
        self.mapper.ref_fastx = get_test_data("samsum_test_2.fasta")
        self.mapper.output_dir = "tests/tmp_mappy_outputs"

        self.max_threads = 4
        self.mapper.threads = self.max_threads

        return

    def tearDown(self) -> None:
        from shutil import rmtree
        if os.path.isdir(self.mapper.output_dir):
            rmtree(self.mapper.output_dir)
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

    def test_mappy_align_multi(self):
        self.mapper.query_fastx = self.pe_query_fastx
        self.mapper.fq_fmt = 0
        ref_hits = self.mapper.mappy_align_multi()
        self.assertEqual(61, len(ref_hits))
        self.assertEqual(1994, sum([len(x) for x in ref_hits.values()]))
        return

    def test_idx_mappy_align(self):
        self.mapper.query_fastx = self.pe_query_fastx
        self.mapper.fq_fmt = 0
        ref_hits = self.mapper.idx_mappy_align()
        self.assertEqual(61, len(ref_hits))
        self.assertEqual(1994, sum([len(x) for x in ref_hits.values()]))
        return


if __name__ == '__main__':
    unittest.main()
