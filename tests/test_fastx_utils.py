import os
import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_fa = get_test_data("fasta_test.fa")
        self.output_dir = "tests/tmp_fastx_utils_tests"

    def tearDown(self) -> None:
        from shutil import rmtree
        if os.path.isdir(self.output_dir):
            rmtree(self.output_dir)
        return

    def test_split_fastx(self):
        from samsum import fastx_utils
        n_parts = 2
        fastx_parts = fastx_utils.split_fastx(self.test_fa,
                                              self.output_dir,
                                              n_parts)
        self.assertEqual(n_parts, len(fastx_parts))
        for f_path in fastx_parts:
            self.assertTrue(os.path.isfile(f_path))
            self.assertTrue(os.stat(f_path).st_size > 10)
        return

    def test_fasta_reader(self) -> None:
        from samsum import fastx_utils
        ref_seq_lengths, _ = fastx_utils.fasta_seq_lengths(fasta_file=self.test_fa)
        self.assertEqual({"contig 1": 16, "contig 2": 25}, ref_seq_lengths)
        return


if __name__ == '__main__':
    unittest.main()
