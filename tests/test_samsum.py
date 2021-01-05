#!/usr/bin/env python

import os
import unittest
import pytest

from .testing_utils import get_test_data


class SamsumTester(unittest.TestCase):
    def setUp(self) -> None:
        self.test_fasta = get_test_data("samsum_test_2.fasta")
        self.test_sam = get_test_data("samsum_test_2.sam")
        self.output_tbl = os.path.join("tests/tmp_table.tsv")
        return

    def tearDown(self) -> None:
        if os.path.isfile(self.output_tbl):
            os.remove(self.output_tbl)
        return

    def test_main(self):
        """ Test whether samsum info runs at all and exits with the write return code """
        from samsum import __main__
        retcode = __main__.main(["samsum"])
        self.assertEqual(1, retcode)
        return

    def test_samsum_info(self):
        from samsum import commands
        retcode = commands.info(["-v"])
        self.assertEqual(0, retcode)
        return

    def test_samsums_stats(self):
        """ Integrative test for samsum stats """
        from samsum import commands
        # Ensure it quits properly with help call
        with pytest.raises(SystemExit):
            commands.stats(["-h"])

        # Test with a normal dataset
        retcode = commands.stats(["--ref_fasta", self.test_fasta,
                                  "--alignments", self.test_sam,
                                  "--output_table", self.output_tbl,
                                  "--aln_percent", str(50),
                                  "--seq_coverage", str(20),
                                  "--map_quality", str(1),
                                  "--sep", "\t"])
        self.assertEqual(0, retcode)
        return


if __name__ == "__main__":
    unittest.main()
