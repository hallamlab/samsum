import unittest
import os
import glob

from samsum import testing_utils as utils


class SamModuleTester(unittest.TestCase):
    def setUp(self) -> None:
        from samsum import classy
        self.unmapped_dat_example = classy.AlignmentDat("UNMAPPED", ['NA', '-530976544', '', '0', '415858.0'])
        self.alignment_dat_example = classy.AlignmentDat("refseq_name", ["read_1", "1", "5S45M", "0", "1.0"])
        self.sam_list = sorted(glob.glob(os.path.join(utils.get_test_data("test-data"), "*.sam")))
        self.refseq = classy.RefSequence("NODE_1", 200)
        sam_dat = {"NODE_1": [["q1", "1", "5S45M", "0", "1.0"],
                              ["q2", "1", "45M", "0", "1.0"],
                              ["q3", "1", "5S145M", "0", "1.0"],
                              ["q4", "1", "5S145M", "0", "1.0"]]}
        for r, alignments in sam_dat.items():
            while alignments:  # type: list
                aln_dat = classy.AlignmentDat(r, alignments.pop())
                self.refseq.alignments.append(aln_dat)
                self.refseq.reads_mapped += 1
        return

    def test_get_mapped_reads(self):
        from samsum import _sam_module
        test_sam = utils.get_test_data('pytest_1.sam')
        self.assertIsInstance(test_sam, str)
        mapping_list = _sam_module.get_mapped_reads(test_sam, False, 10, 0, 'q')
        self.assertEqual(8, len(mapping_list))
        return

    def test_load_sam(self):
        test_aln_data = ["query_read_name", "1", "5S145M", "0", "1.0"]
        self.alignment_dat_example.load_sam(test_aln_data)
        self.assertEqual(1, self.alignment_dat_example.start)
        self.assertEqual(1.0, self.alignment_dat_example.weight)
        self.assertEqual("refseq_name", self.alignment_dat_example.ref)
        self.assertEqual(145, self.alignment_dat_example.end)
        self.assertEqual(145, self.alignment_dat_example.decode_cigar())
        return

    def test_load_unmapped(self):
        self.assertEqual(415858.0, self.unmapped_dat_example.weight)
        self.assertEqual("UNMAPPED", self.unmapped_dat_example.ref)
        return

    def test_decode_cigar(self):
        cigar_str_1 = "101S19M30S"
        self.alignment_dat_example.cigar = cigar_str_1
        self.assertEqual(19, self.alignment_dat_example.decode_cigar())
        self.assertEqual(150, self.alignment_dat_example.read_length)
        cigar_str_2 = "101S19M30H"
        self.alignment_dat_example.cigar = cigar_str_2
        self.assertEqual(19, self.alignment_dat_example.decode_cigar())
        self.assertEqual(120, self.alignment_dat_example.read_length)

        return

    def test_ref_sequence_abundances(self):
        from samsum import commands
        from samsum.classy import RefSequence
        test_sam = utils.get_test_data("samsum_test_2.sam")
        test_asm = utils.get_test_data("samsum_test_2.fasta")
        ref_seq_abunds = commands.ref_sequence_abundances(aln_file=test_sam, seq_file=test_asm,
                                                          min_aln=10, p_cov=0, map_qual=0, multireads=True)
        e10 = ref_seq_abunds['AB-755_P22_E10_NODE_6_length_36342_cov_2156.57_ID_21']  # type: RefSequence
        self.assertEqual(e10.reads_mapped, 10)
        self.assertEqual(e10.weight_total, 5.0)

        # Test the number of reads that mapped by summing their total weights
        total_mapped = 0.0
        for refseq in ref_seq_abunds.values():  # type: RefSequence
            total_mapped += refseq.reads_mapped
        self.assertEqual(220.0, total_mapped)

        self.assertTrue(1E6-1 < sum(refseq.tpm for refseq in ref_seq_abunds.values()) < 1E6+1)
        return

    def test_proportion_covered(self):
        self.assertEqual(4, self.refseq.reads_mapped)
        self.assertEqual(0.72, self.refseq.proportion_covered())
        return


if __name__ == "__main__":
    unittest.main()
